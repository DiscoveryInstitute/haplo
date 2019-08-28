#include "FileSNPsWriter.h"
#include <algorithm>
#include <fstream>
#include <iostream>
using namespace std;

void FileSNPsWriter::WriteSNPs(
        const string& output_filename,
        const HaploBlockVariants& variants,
        const pair<vector<PersonNode>, vector<PersonNode>>& whole_generation,
        mt19937_64& rng) {

    const long total_blocks = variants.getNumberBlocks();
    const auto& boundary_loci = variants.getBlockBoundaries();
    const long total_males = whole_generation.first.size();
    const long total_females = whole_generation.second.size();
    const long total_chromosomes = 2 * (total_males + total_females);

    ofstream ofs(output_filename, ofstream::out | ofstream::app);
    ofs << "NAMES";
    for (int i=0; i<total_chromosomes; ++i) ofs << '\t' << i;
    ofs << endl;
    ofs << "REGION\tchr\t0\t" + to_string(boundary_loci.back()) << endl;

    // For every haplotype block.
    for (long ib=0; ib<total_blocks; ++ib) {
        LocusId block_locus = boundary_loci[ib];
        LocusId last_locus = boundary_loci[ib+1]-1;
        const MutationGraph& graph = variants.getMutationGraph(block_locus);

        // Get the position/allele pairs
        const auto allele_ids = graph.getAlleleIds();
        vector<pair<LocusId,AlleleId>> mutation_positions_alleles;
        for (AlleleId allele_id : allele_ids) {
            const auto nmuts = graph.getNumberOfMutations(allele_id);
            const auto& muts = graph.getMutations(allele_id);
            if ((size_t)nmuts == muts.size()) {
                for (int i=0; i<nmuts; ++i) {
                    LocusId mut = muts[i].first;
                    mutation_positions_alleles.emplace_back(mut,allele_id);
                }
            }
            else {
                for (int i=0; i<nmuts; ++i) {
                    LocusId mut = uniform_int_distribution<long>(block_locus,last_locus)(rng);
                    mutation_positions_alleles.emplace_back(mut, allele_id);
                }
            }
        }
        // Sort and space them, throw away any that don't fit.
        auto& vec = mutation_positions_alleles;
        sort(vec.begin(), vec.end());
        const int nsites = vec.size();
        const long block_length = last_locus - block_locus + 1;
        if ((size_t)block_length<vec.size()) vec.resize(block_length);
        for (int i=1; i<nsites; ++i) {
            if (vec[i-1].first>=vec[i].first) {
                vec[i].first=vec[i-1].first + 1;
            }
            while (vec[i].first>last_locus) {
                vec[i].first--;
                for (int j=i-1; j>=0; --j) {
                    if (vec[j].first==vec[j+1].first) vec[j].first--;
                }
            }
        }

        // For every position in order on the block.
        for (auto& entry : mutation_positions_alleles) {
            const LocusId position = entry.first;
            const AlleleId source_allele_id = entry.second;

            string line(total_chromosomes, ' ');
            int ic=0;
            for (auto* persons : { &whole_generation.first, &whole_generation.second })
            for (auto& person : *persons)
            for (auto* chromosome : { &person.chromosomeFromFather, &person.chromosomeFromMother }) {
                AlleleId allele = chromosome->getValue(block_locus);
                bool variant_present = graph.isDerivativeAllele(source_allele_id, allele);
                line[ic++] = variant_present ? 'A' : 'C';
            }
            ofs << position << '\t' << line << endl;
        }
    }

}

