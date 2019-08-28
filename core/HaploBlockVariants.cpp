#include <algorithm>
#include <iterator>
#include <random>
#include <vector>
#include "Chromosome.h"
#include "HaploBlockBoundaries.h"
#include "HaploBlockVariants.h"

using namespace std;

vector<LocusId> generateMutations(int nmuts, LocusId start_locus, LocusId end_locus, mt19937_64& rng) {
    uniform_int_distribution<long> distribution(start_locus, end_locus-1);
    vector<LocusId> mutations; mutations.reserve(nmuts);
    for (long i=0; i<nmuts; i++) mutations.push_back(distribution(rng));
    sort(mutations.begin(), mutations.end()); // sort
    return mutations;
}


/**
 *  Each original allele is represented by a set of single nucleotide changes from some
 *  (unspecified) reference sequence.
 *  Each original allele may be placed on more than one original chromosome.
 * */
struct OriginalAllele {
    vector<LocusId> single_nucleotide_differences;
    vector<int> chromosome_placements;
};

HaploBlockVariants::HaploBlockVariants(
        const HaploBlockBoundaries& hbb,
        const long max_nodes_per_block,
        const long primordial_block_length,
        const double probability_dimorphic,
        const double probability_tetramorphic,
        const double primordial_diversity,
        const bool use_mutation_loci,
        mt19937_64& rng)
{
    m_max_nodes_per_block = max_nodes_per_block;
    LocusId chromosome_length = hbb.getChromosomeLength();
    long n_primordial_blocks =
            chromosome_length / primordial_block_length
            + (chromosome_length % primordial_block_length ? 1 : 0);
    double primordial_snp_rate =  // primordial mutation-like SNPs
            primordial_diversity /
            (2 * probability_dimorphic + 4 * probability_tetramorphic);

    vector<LocusId> primordial_block_boundaries = hbb.getRandomBlockBoundaries(n_primordial_blocks, rng);
    vector<vector<OriginalAllele>> primordial_blocks(n_primordial_blocks);
    for (long ipb=0; ipb<n_primordial_blocks; ++ipb) {
        LocusId block_start = primordial_block_boundaries.at(ipb);
        LocusId block_end = primordial_block_boundaries.at(ipb+1);
        double r = uniform_real_distribution<double>(0,1.0)(rng);
        int n_original_alleles = (r-=probability_tetramorphic) < 0 ? 4 :
                                 (r-=probability_dimorphic) < 0  ? 2 : 1;
        int n_copies_per_original_allele = 4 / n_original_alleles;
        vector<int> placements{0,1,2,3};
        swap(placements.at(0), placements.at(uniform_int_distribution<int>(0,3)(rng)));
        swap(placements.at(1), placements.at(uniform_int_distribution<int>(1,3)(rng)));
        swap(placements.at(2), placements.at(uniform_int_distribution<int>(2,3)(rng)));
        auto it = placements.begin();
        for (int ioa=0; ioa<n_original_alleles; ++ioa) {
            vector<LocusId> snps;
            vector<int> plcs;
            if (n_original_alleles > 1) {
                double expected_n_snps = primordial_snp_rate * (block_end - block_start);
                long n_snps = poisson_distribution<long>(expected_n_snps)(rng);
                snps = generateMutations(n_snps, block_start, block_end, rng);
            }
            for (int i=0; i<n_copies_per_original_allele; ++i) plcs.push_back(*it++);
            OriginalAllele original_allele { snps, plcs };
            primordial_blocks.at(ipb).push_back(original_allele);
        }

    }

    const auto& recombination_block_boundaries = hbb.getRecordedBlockBoundaries();
    set<LocusId> all_blocks;
    for (LocusId boundary : primordial_block_boundaries) all_blocks.insert(boundary);
    for (LocusId boundary : recombination_block_boundaries) all_blocks.insert(boundary);
    m_boundaries = vector<LocusId>(all_blocks.begin(), all_blocks.end());

    long n_all_blocks = m_boundaries.size() - 1;
    m_variants.reserve(n_all_blocks);

    vector<vector<pair<LocusId,AlleleId>>> reference_chromosomes(4);

    for (long ib=0; ib<n_all_blocks; ++ib) {

        LocusId block_start = m_boundaries.at(ib);
        LocusId block_end = m_boundaries.at(ib+1);
        const auto& pbb = primordial_block_boundaries;
        long ipb = distance(pbb.begin(), upper_bound(pbb.begin(), pbb.end(), block_start)) - 1L; // index of primordial block
        vector<vector<LocusId>> observed_block_founding_snd_sets;
        AlleleId allele_id = 1;
            for (const auto& original_allele : primordial_blocks.at(ipb)) {
            vector<LocusId> snd_set;
            for (LocusId snd : original_allele.single_nucleotide_differences) {
                if (block_start <= snd && snd < block_end) snd_set.push_back(snd);
            }
            observed_block_founding_snd_sets.push_back(snd_set);
            for (int iref : original_allele.chromosome_placements) {
                reference_chromosomes.at(iref).emplace_back(block_start, allele_id);
            }
            ++allele_id;
        }
        m_variants.emplace_back(observed_block_founding_snd_sets);
        m_block_ids.emplace(block_start, ib);  // maps locus_id -> block_id

    }
    m_block_ids.emplace(chromosome_length, n_all_blocks);

    for (auto& reference_chromosome : reference_chromosomes) {
        Chromosome new_chromosome;
        new_chromosome.m_intervals = move(reference_chromosome);
        new_chromosome = new_chromosome.compressed();
        m_reference_chromosomes.push_back(new_chromosome);
    }

}


const Chromosome& HaploBlockVariants::getReferenceChromosome(long i) {
    long n = m_reference_chromosomes.size();
    return m_reference_chromosomes.at(i % n);
}


Chromosome HaploBlockVariants::mutate(
        const Chromosome &chromosome_in, double nucleotide_mutation_rate, long generation, mt19937_64 &rng) {

    LocusId chromosome_length = m_boundaries.back();
    double total_mutation_rate = nucleotide_mutation_rate * chromosome_length;
    int n_mutations = poisson_distribution<int>(total_mutation_rate)(rng);
    if (n_mutations == 0) return chromosome_in;

    const long end_locus = m_boundaries.back();
    vector<LocusId> mutationLoci = generateMutations(n_mutations, 0, end_locus, rng);

    map<BlockId, vector<LocusId>> mutatedBlocks;
    for (LocusId locus :  mutationLoci) {
        BlockId block = prev(m_block_ids.upper_bound(locus))->second;
        mutatedBlocks[block].push_back(locus);
    }

    Chromosome chromosome_out;
    auto& out = chromosome_out.m_intervals;

    auto begin_it = chromosome_in.m_intervals.cbegin();
    auto end_it = chromosome_in.m_intervals.cend();
    auto it = begin_it;
    AlleleId allele = UNKNOWN_ALLELE;

    for (auto block_it = mutatedBlocks.begin(); block_it != mutatedBlocks.end(); ++block_it) {
        BlockId block = block_it->first;
        auto& muts = block_it->second;
        LocusId block_begin = m_boundaries.at(block);
        LocusId block_end = m_boundaries.at(block+1);

        // Include interval boundaries that occur before the start of the mutation block.
        for ( ; it != end_it && it->first <= block_begin; ++it) {
            if (it->second != allele) out.push_back(*it);
            allele = it->second;
        }
        if (allele == UNKNOWN_ALLELE) continue; // Ignore unknown/non-ancestral

        // Create a mutated allele.
        AlleleId mutated_allele =
                m_variants.at(block).addAlleleMutation(allele, generation, muts);

        // Add the boundary for the start of the mutated block
        //   (or if there is an existing boundary, update it, or even remove it).
        if (out.back().first < block_begin) {
            out.emplace_back(block_begin, mutated_allele);   // add entry

        } else {
            out.back().second = mutated_allele; // update last entry
            if (out.size() >= 2 && out.end()[-2].second == mutated_allele) {
                out.pop_back();                 // remove last entry
            }
        }

        // Add a boundary for the end of the mutated block, if necessary.
        auto interval_end = it != end_it ? it->first : end_locus;
        if (block_end < interval_end) {
            out.emplace_back(block_end, allele);
        } else {
            allele = mutated_allele;
        }

    }
    // Include interval boundaries that occur after the last mutated block.
    for ( ; it != end_it; ++it) {
        if (it->second != allele) out.push_back(*it);
        allele = it->second;
    }

    return chromosome_out;

}


bool HaploBlockVariants::sampleChromosomesAndSimplify(
        const pair<vector<PersonNode>, vector<PersonNode>>& whole_generation, bool hard) {

    bool simplified = false;
    const long end_locus = m_boundaries.back();
    for (LocusId locus : m_boundaries) {
        if (locus == end_locus) break;
        auto bi = m_block_ids.at(locus);
        MutationGraph& mut_graph = m_variants.at(bi);
        if (!hard && !mut_graph.tooManyNodes(m_max_nodes_per_block)) continue;
        mut_graph.clearAlleleFrequencies();
        for (auto* persons : { &whole_generation.first, &whole_generation.second })
        for (auto& person : *persons)
        for (auto* chromosome : { &person.chromosomeFromFather, &person.chromosomeFromMother }) {
            AlleleId allele = chromosome->getValue(locus);
            if (allele == UNKNOWN_ALLELE) continue;
            mut_graph.incrementAlleleFrequency(allele);
        }
        mut_graph.simplify();
        simplified = true;
    }
    return simplified;

}

Chromosome HaploBlockVariants::getFounderAlleleIdChromosome(
        const Chromosome& normal_id_chromosome) const {
    Chromosome new_chromosome;
    new_chromosome.m_intervals = normal_id_chromosome.m_intervals;
    for (auto& interval : new_chromosome.m_intervals) {
        LocusId locus = interval.first;
        AlleleId& allele = interval.second;
        if (allele==UNKNOWN_ALLELE) continue;
        BlockId block = prev(m_block_ids.upper_bound(locus))->second;
        AlleleId founder_allele = m_variants.at(block).getFounderId(allele);
        allele = (allele==founder_allele) ? founder_allele : -founder_allele;
    }
    return new_chromosome;
}

const MutationGraph& HaploBlockVariants::getMutationGraph(LocusId locus) const {
    return m_variants.at(m_block_ids.at(locus));
}

void HaploBlockVariants::deleteMutationGraph(LocusId locus) {
    m_variants.at(m_block_ids.at(locus)) = MutationGraph();
}

bool HaploBlockVariants::operator==(const HaploBlockVariants& o) const {
    return m_max_nodes_per_block   == o.m_max_nodes_per_block
        && m_reference_chromosomes == o.m_reference_chromosomes
        && m_variants              == o.m_variants
        && m_boundaries            == o.m_boundaries
        && m_block_ids             == o.m_block_ids;
}
