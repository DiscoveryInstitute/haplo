#pragma once
#include <map>
#include <random>
#include <utility>
#include <vector>
#include "Chromosome.h"
#include "HaploBlockBoundaries.h"
#include "MutationGraph.h"
#include "PersonNode.h"


/**
    HaploBlockVariants:
    Stores all the alleles / haplotype variants that exist on each haplotype block
    (non-recombining section of chromosome), and any tree-like relationship between them.
    Generates new alleles by mutation, recording which allele it mutated from.
    Also stores four reference chromosomes (the chromosomes of the founding pair).
    Also stores boundaries between blocks - copy of information in HaploBlockBoundaries.
*/
class HaploBlockVariants {

        friend class FileCheckpointer;

    private:

        long m_max_nodes_per_block;
        std::vector<Chromosome> m_reference_chromosomes;
        std::vector<MutationGraph> m_variants;
        std::vector<LocusId> m_boundaries;         // map BlockId -> LocusId
        std::map<LocusId, BlockId> m_block_ids;    // map LocusId -> BlockId

    public:

        HaploBlockVariants() : m_max_nodes_per_block(0L) {} // for creating unassigned variables
        bool unassigned() const { return m_reference_chromosomes.empty() && m_variants.empty() && m_boundaries.empty(); }

        HaploBlockVariants(const HaploBlockBoundaries& hbb,
                           const long max_nodes_per_block,
                           const long primordial_block_length,
                           const double probability_dimorphic,
                           const double probability_tetramorphic,
                           const double primordial_diversity,
                           const bool use_mutation_loci,
                           std::mt19937_64& rng);

        LocusId getChromosomeLength() const { return m_boundaries.back(); }
        long getNumberBlocks() const { return m_boundaries.size()-1; }
        const std::vector<LocusId> getBlockBoundaries() const { return m_boundaries; }

        const Chromosome& getReferenceChromosome(long chromosomeIndex);

        Chromosome mutate(const Chromosome&, double mu, long gen, std::mt19937_64& rng);

        bool sampleChromosomesAndSimplify(
                const std::pair<std::vector<PersonNode>, std::vector<PersonNode>>& whole_gen, bool hard = false);

        Chromosome getFounderAlleleIdChromosome(const Chromosome& normal_id_chromosome) const;

        const MutationGraph& getMutationGraph(LocusId locus) const;
        void deleteMutationGraph(LocusId locus);

        bool operator==(const HaploBlockVariants&) const;
};
