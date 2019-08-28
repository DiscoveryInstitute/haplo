#include "Analysis.h"
#include <map>
#include <utility>
#include <vector>
#include "Chromosome.h"
#include "GraphDataBinned.h"
#include "GraphDataBinned2D.h"
#include "Parser.h"

using namespace std;

#pragma omp declare reduction (merge : GraphDataBinned : omp_out.merge(omp_in) ) initializer (omp_priv = omp_orig.zeroedCopy())
#pragma omp declare reduction (merge2D : GraphDataBinned2D : omp_out.merge(omp_in) ) initializer (omp_priv = omp_orig.zeroedCopy())

/**
 *  Calculate the output parameters listed in table 4, page 14 of
 *     http://bio-complexity.org/ojs/index.php/main/article/view/BIO-C.2016.4
 *  The method assumes that all mutation_graphs have already been pruned.
 */
Analysis Analysis::Calculate(
        Logger& logger,
        Parser& inputs,
        HaploBlockVariants& variants,
        pair<vector<PersonNode>, vector<PersonNode>>& whole_generation) {

    logger.logEvent("Analysis begins.");
    const long max_bins = inputs.getPositiveLong("analysis_datapoints_maximum");
    const long max_distance = inputs.getPositiveLong("analysis_linkage_distance_maximum");
    const bool do_linkage_stats = inputs.getBool("analysis_do_linkage_stats");
    const bool use_mutation_loci = inputs.getBool("use_mutation_loci");
    const double min_freq = inputs.getNonNegativeDouble("analysis_linkage_minimum_frequency");

    const long total_sn_loci = variants.getChromosomeLength();
    const long total_blocks = variants.getNumberBlocks();
    const auto& boundary_loci = variants.getBlockBoundaries();
    const LocusId chromosome_length = boundary_loci.back();
    const auto block_loci = vector<LocusId>(boundary_loci.begin(),boundary_loci.end()-1);
    const long total_males = whole_generation.first.size();
    const long total_females = whole_generation.second.size();
    const long total_chromosomes = 2 * (total_males + total_females);
    const AlleleId NO_PARENT = -1;
    const Chromosome EMPTY = Chromosome();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Sample and Simplify Mutation Graphs
    variants.sampleChromosomesAndSimplify(whole_generation, true);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Transform the data into a more cache-friendly shape to make processing faster.
    logger.logEvent(" Data Transform begins.");
    vector<vector<AlleleId>> allele_data(total_blocks);      // [BlockId][ChromosomeId]
    vector<vector<AlleleId>> parent_data(total_blocks);      // [BlockId][AlleleId]
    vector<vector<long>> nsnp_data(total_blocks);             // [BlockId][AlleleId];
    vector<vector<vector<LocusId>>> snp_data;                // [BlockId][AlleleId][int]
    if (use_mutation_loci) snp_data.resize(total_blocks);

    // Extract and transpose the chromosome data.
    bool empty_chromosomes = false;
    for (BlockId ib=0; ib<total_blocks; ++ib) allele_data[ib].reserve(total_chromosomes);
    for (auto* persons : { &whole_generation.first, &whole_generation.second })
    for (auto& person : *persons)
    for (auto* chromosome : { &person.chromosomeFromFather, &person.chromosomeFromMother }) {
        if (*chromosome == EMPTY) { empty_chromosomes = true; continue; }
        #pragma omp parallel for
        for (BlockId bi=0; bi<total_blocks; ++bi) {
            LocusId block_locus = block_loci[bi];
            AlleleId allele = chromosome->getValue(block_locus);
            if (allele==UNKNOWN_ALLELE) throw logic_error("Unknown allele in non-empty chromosome.");
            allele_data[bi].push_back(allele);
        }
        *chromosome = EMPTY; // delete/free memory
    }
    whole_generation.first = {}; // delete/free memory
    whole_generation.second = {}; // delete/free memory
    if (empty_chromosomes) logger.logEvent("  There were empty chromosomes in sample.");
    logger.logEvent("  Data Transposed.");

    // Next, map the non-contiguous allele ids in MutationGraph to a set of contiguous allele ids.
    #pragma omp parallel for
    for (BlockId bi=0; bi<total_blocks; ++bi) {
        const MutationGraph& graph = variants.getMutationGraph(block_loci[bi]);
        // Get old ids. Map to new ids.
        vector<AlleleId> old_allele_ids = graph.getAlleleIds();
        map<AlleleId, AlleleId> new_allele_id_map;
        new_allele_id_map[UNKNOWN_ALLELE] = NO_PARENT;
        AlleleId new_allele_id = 0;
        for (const auto& old_allele_id : old_allele_ids) {
            new_allele_id_map[old_allele_id] = new_allele_id++;
        }
        // Copy the (old) parent_ids and mutation-locus groups into a contiguous vector structure.
        int n_alleles = old_allele_ids.size();
        parent_data[bi] = vector<AlleleId>(n_alleles);
        nsnp_data[bi] = vector<long>(n_alleles);
        if (use_mutation_loci) snp_data[bi] = vector<vector<LocusId>>(n_alleles);

        for (const auto& old_allele_id : old_allele_ids) {
            AlleleId new_allele_id = new_allele_id_map.at(old_allele_id);
            parent_data[bi][new_allele_id] = graph.getParentId(old_allele_id);
            // important that the graph was pruned first ^

            nsnp_data[bi][new_allele_id] = graph.getNumberOfMutations(old_allele_id);
            if (use_mutation_loci) {
                const auto& mutations = graph.getMutations(old_allele_id);
                snp_data[bi][new_allele_id].reserve(mutations.size());
                for (auto& entry : mutations) snp_data[bi][new_allele_id].push_back(entry.first);
            }
        }

        // Finally update all the copied allele ids from old to new.
        for (AlleleId& allele_id : allele_data[bi]) allele_id = new_allele_id_map.at(allele_id);
        for (AlleleId& parent_id : parent_data[bi]) parent_id = new_allele_id_map.at(parent_id);

        variants.deleteMutationGraph(block_loci[bi]); // delete/free memory
    }
    variants = HaploBlockVariants(); // delete/free memory
    logger.logEvent(" Data Transform complete.");


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Process the data to get statistics.

    const double total_freq = (double)allele_data[0].size();
    logger.logEvent("Total samples: " + to_string(total_freq));

    const double freq_max = 0.5;
    const size_t freq_nbins = min(max_bins, (long)ceil(freq_max * total_freq));
    const double link_max = min(total_sn_loci, max_distance);
    const size_t link_nbins = min(max_bins, max_distance);

    long total_snp_sites = 0L;
    double sn_diversity = 0.0;
    GraphDataBinned snp_freq_spectrum(0, freq_max, freq_nbins, true);
    GraphDataBinned snp_distribution(0, chromosome_length, min(chromosome_length,max_bins));
    GraphDataBinned correlation_function(0, link_max, link_nbins);
    GraphDataBinned dprime_function(0, link_max, link_nbins);
    GraphDataBinned complete_linkage_function(0, link_max, link_nbins);
    GraphDataBinned almost_complete_linkage_function(0, link_max, link_nbins);
    GraphDataBinned normalized_kullback_leibler_function(0, link_max, link_nbins);
    GraphDataBinned sigma_squared_numerator(0, link_max, link_nbins);
    GraphDataBinned sigma_squared_denominator(0, link_max, link_nbins);
    GraphDataBinned2D freq_pair_spectrum({0.0,0.5},{0.0,0.5}, 20, true);

    // First, the single-locus statistics.
    logger.logEvent(" Single Nucleotides statistics begin.");
    long n_alleles_fixed=0;
    #pragma omp parallel for reduction(+: total_snp_sites, sn_diversity) \
                             reduction(merge: snp_freq_spectrum)
    for (BlockId bi=0; bi<total_blocks; ++bi) {

        const int n_alleles = nsnp_data[bi].size();
        double freqs[n_alleles];                    // Roughly corresponds to pi in the paper.
        for (AlleleId ai=0; ai<n_alleles; ++ai) freqs[ai] = 0.0;

        // First, the frequencies for all alleles (unique combinations of single nucleotides).
        for (AlleleId ai : allele_data[bi]) {
            ++freqs[ai];
        }

        // Second, transform these into frequencies for the single nucleotides.
        for (AlleleId ai=n_alleles-1; ai>=0; --ai) {
            AlleleId pai = parent_data[bi][ai];
            if (pai == NO_PARENT) continue;
            freqs[pai] += freqs[ai];
        }

        // Finally, update the single nucleotide stats.
        double norm_sn_diversity = 2.0 / total_sn_loci;
        for (AlleleId ai=0; ai<n_alleles; ++ai) {
            double fa = freqs[ai] / total_freq;
            if (fa==0.0) throw logic_error("freqs elements should not be empty");
            if (fa>=1.0) { ++n_alleles_fixed; continue; }
            fa = fa < 0.5 ? fa : (1.0-fa);
            long n_snps = nsnp_data[bi][ai];
            total_snp_sites += n_snps;
            sn_diversity += norm_sn_diversity * n_snps * fa * (1.0 - fa);
            snp_freq_spectrum.sampleWeight(fa, n_snps);
            if (use_mutation_loci) {
                for (auto snp : snp_data[bi][ai]) snp_distribution.sampleWeight(snp,1);
            }
        }
    }
    logger.logEvent("  Fixed Alleles: " + to_string(n_alleles_fixed));
    logger.logEvent(" Single Nucleotides statistics complete.");

    if (do_linkage_stats) {

    // Second, the two-locus statistics.
    const size_t nblocks = block_loci.size();
    const double f1 = 1.0 / total_freq;
    logger.logEvent(" Pair Nucleotide statistics begin.");
    logger.logEvent("  Total blocks: " + to_string(nblocks) + ".");
    #pragma omp parallel for schedule(dynamic) \
                             reduction(merge: correlation_function, \
                                              dprime_function, \
                                              complete_linkage_function, \
                                              almost_complete_linkage_function, \
                                              normalized_kullback_leibler_function, \
                                              sigma_squared_numerator, \
                                              sigma_squared_denominator) \
                             reduction(merge2D: freq_pair_spectrum)
    for (size_t bi=0; bi<nblocks; ++bi) {
        for (size_t bj=bi; bj<nblocks; ++bj) {

            const LocusId block_locus_i = block_loci[bi];
            const LocusId block_locus_j = block_loci[bj];
            const LocusId block_end_locus_i =
                    (bi+1)==nblocks ? total_sn_loci : block_loci[bi+1];
            const long block_distance = block_locus_j - block_end_locus_i;
            if (block_distance > max_distance) break; // goto next bi

            const int n_alleles_i = nsnp_data[bi].size();
            const int n_alleles_j = nsnp_data[bj].size();
            double freqs_i[n_alleles_i];
            double freqs_j[n_alleles_j];
            double freqs_ij[n_alleles_i][n_alleles_j];
            for (AlleleId ai=0; ai<n_alleles_i; ++ai) freqs_i[ai] = 0.0;
            for (AlleleId aj=0; aj<n_alleles_j; ++aj) freqs_j[aj] = 0.0;
            for (AlleleId ai=0; ai<n_alleles_i; ++ai)
            for (AlleleId aj=0; aj<n_alleles_j; ++aj) freqs_ij[ai][aj] = 0.0;

            // First, the frequencies for all pairs of alleles in block i and block j.
            for (ChromosomeId ci=0; ci<total_freq; ++ci) {
                AlleleId ai = allele_data[bi][ci];
                AlleleId aj = allele_data[bj][ci];
                ++freqs_i[ai];
                ++freqs_j[aj];
                ++freqs_ij[ai][aj];
            }
            // Second, transform these into frequencies for the single nucleotides, for block i ...
            for (AlleleId ai=n_alleles_i-1; ai>=0; --ai) {
                AlleleId pai = parent_data[bi][ai];
                if (pai == NO_PARENT) continue;
                freqs_i[pai] += freqs_i[ai];
                for (AlleleId aj=0; aj<n_alleles_j; ++aj) {
                    freqs_ij[pai][aj] += freqs_ij[ai][aj];
                }
            }
            // ... then for block j
            for (AlleleId aj=n_alleles_j-1; aj>=0; --aj) {
                AlleleId paj = parent_data[bj][aj];
                if (paj == NO_PARENT) continue;
                freqs_j[paj] += freqs_j[aj];
                for (AlleleId ai=0; ai<n_alleles_i; ++ai) {
                    freqs_ij[ai][paj] += freqs_ij[ai][aj];
                }
            }

            // Calculate correlations for this pair of blocks
            double kl_numerator = 0.0;
            double kl_denom_sqrd = 0.0;
            for (AlleleId ai=0; ai<n_alleles_i; ++ai) {
                double fa = freqs_i[ai] / total_freq;
                if (fa==0.0) throw logic_error("freqs elements should not be empty");
                if (fa==1.0) continue;
                if (fa<min_freq || fa>1.0-min_freq) continue;

                for (AlleleId aj=(bi==bj?ai:0); aj<n_alleles_j; ++aj) {
                    double fb = freqs_j[aj] / total_freq;
                    if (fb==0.0) throw logic_error("freqs elements should not be empty");
                    if (fb==1.0) continue;
                    if (fb<min_freq || fb>1.0-min_freq) continue;

                    double fab = freqs_ij[ai][aj] / total_freq;
                    double fA = 1 - fa;
                    double fB = 1 - fb;
                    double faB = fa - fab;
                    double fAb = fb - fab;
                    double fAB = 1 - fab - faB - fAb;
                    double delta = fab * fAB - faB * fAb;
                    double sig2_numerator = delta * delta;
                    double sig2_denominator = (fa * fA * fb * fB);
                    double correlation = sig2_numerator / sig2_denominator; // eqn 47a: r^2
                    double dprime = delta / (delta > 0 ? min(fa*fB,fA*fb) : -min(fa*fb,fA*fB));
                    double complete_linkage        = (fab==0  || faB==0  || fAb==0  || fAB == 0) ? 1.0 : 0.0;
                    double almost_complete_linkage = (fab<=f1 || faB<=f1 || fAb<=f1 || fAB <=f1) ? 1.0 : 0.0;

                    double mafA = min(fa,fA);
                    double mafB = min(fb,fB);

                    if (use_mutation_loci) {
                        for (LocusId locus_i : snp_data[bi][ai])
                        for (LocusId locus_j : snp_data[bj][aj]) {
                            double locus_distance = fabs(locus_i - locus_j);
                            correlation_function.sampleValue(locus_distance, correlation);
                            dprime_function.sampleValue(locus_distance, dprime);
                            complete_linkage_function.sampleValue(locus_distance, complete_linkage);
                            almost_complete_linkage_function.sampleValue(locus_distance, almost_complete_linkage);
                            sigma_squared_numerator.sampleValue(locus_distance, sig2_numerator);
                            sigma_squared_denominator.sampleValue(locus_distance, sig2_denominator);
                            if (0<locus_distance && locus_distance < 1e4) {
                                freq_pair_spectrum.sampleWeight(mafA, mafB);
                                freq_pair_spectrum.sampleWeight(mafB, mafA);
                            }
                        }
                    } else {
                        double distance = fabs(block_locus_j-block_locus_i);
                        long nlocuspairs = nsnp_data[bi][ai] * nsnp_data[bj][aj];
                        correlation_function.sampleValue(             distance, correlation,             nlocuspairs );
                        dprime_function.sampleValue(                  distance, dprime,                  nlocuspairs );
                        complete_linkage_function.sampleValue(        distance, complete_linkage,        nlocuspairs );
                        almost_complete_linkage_function.sampleValue( distance, almost_complete_linkage, nlocuspairs );
                        sigma_squared_numerator.sampleValue(          distance, sig2_numerator,          nlocuspairs );
                        sigma_squared_denominator.sampleValue(        distance, sig2_denominator,        nlocuspairs );
                        if (0 < distance && distance < 1e4) {
                            freq_pair_spectrum.sampleWeight(mafA, mafB, nlocuspairs);
                            freq_pair_spectrum.sampleWeight(mafB, mafA, nlocuspairs);
                        }
                    }

                    bool ai_is_founder = parent_data[bi][ai] == NO_PARENT;
                    bool aj_is_founder = parent_data[bj][aj] == NO_PARENT;
                    if (ai_is_founder && aj_is_founder) {
                        double weight = (bi==bj && ai!=aj) ? 2.0 : 1.0;
                        kl_numerator += weight * (fab ? fab * log(fab / (fa * fb)) : 0.0);
                        kl_denom_sqrd += weight * (fa * fb * log(fa) * log(fb));
                    }
                }
            }
            if (kl_denom_sqrd) {
                double kl_distance = block_locus_j - block_locus_i;
                double kl_value = kl_numerator / sqrt(kl_denom_sqrd);
                normalized_kullback_leibler_function.sampleValue(kl_distance, kl_value);
            }
            logger.logProgress("  Pair " + to_string(bi) + ", " + to_string(bj) + " done.");
        } // end for block bj
    } // end for block bj
    logger.logEvent(" Pair Nucleotide statistics complete.");

    sigma_squared_numerator.divide(sigma_squared_denominator);

    } // end if do_linkage

    logger.logEvent("Analysis complete.");

    return Analysis {
        total_blocks,
        total_snp_sites,
        sn_diversity,
        move(snp_freq_spectrum),
        move(snp_distribution),
        move(correlation_function),
        move(dprime_function),
        move(complete_linkage_function),
        move(almost_complete_linkage_function),
        move(normalized_kullback_leibler_function),
        move(sigma_squared_numerator),
        move(freq_pair_spectrum)
    };

}

