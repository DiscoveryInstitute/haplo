#pragma once

#include <map>
#include <utility>
#include <vector>
#include "Chromosome.h"

/**
    MutationGraph:
    Stores a set of variants of a haplotype block, each variant referenced by an allele id.
    For each allele, stores the allele from which it was derived, plus the mutations that makes
    it different (or number of mutations if use_mutation_loci is not known), and the number of
    alleles that are further derived from itself, forming a tree graph.
    Can also be used to count the number of copies of each allele in the ancestral population,
    for summarizing statistics, or for culling the tree for efficiency reasons.
*/
class MutationGraph {

        friend class FileCheckpointer;

	private:

        // Represents an allele or a particular subset of mutations that appear together.
        struct Node {
            // The allele from which this one was derived.
            AlleleId parent_id;
            // The mutations (SNPs - Single Nucleotide Polymorphisms), represented by locations
            //   on the chromosomes and age in generations since they occurred.
            long n_mutations;
            std::vector<std::pair<LocusId, long>> mutations;
            // The number of alleles that are (directly) derived from this one.
            //   (and thus also have these SNPs).
            int n_child_nodes;
            // The number of copies of this particular allele / SNP group exist in the population.
            long frequency;

            bool operator==(const Node&) const;
        };

        bool m_use_mutation_loci;
        std::map<AlleleId, Node> m_nodes;
        AlleleId m_next_id;

	public:

        MutationGraph() : m_use_mutation_loci(false), m_next_id(0L) {} // for creating unassigned variables
        MutationGraph(std::vector<std::vector<LocusId>> founder_alleles, const bool use_mutation_loci=true);

        AlleleId addAlleleMutation(AlleleId parent_id, long gen, const std::vector<LocusId>& loci);

        inline bool tooManyNodes(size_t nnodes) const { return m_nodes.size() > nnodes; }
        inline void incrementAlleleFrequency(AlleleId id) { m_nodes.at(id).frequency++; }
        void simplify();
        inline void clearAlleleFrequencies() { for (auto& entry : m_nodes) entry.second.frequency = 0L; }

        std::vector<AlleleId> getAlleleIds() const;
        AlleleId getParentId(AlleleId id) const { return m_nodes.at(id).parent_id; }
        AlleleId getFounderId(AlleleId id) const;
        bool isDerivativeAllele(AlleleId parent_id, AlleleId id) const;
        const std::vector<std::pair<LocusId,long>>& getMutations(AlleleId id) const { return m_nodes.at(id).mutations; }
        int getNumberOfMutations(AlleleId id) const { return m_nodes.at(id).n_mutations; }
        int getMutationDistance(AlleleId ia, AlleleId ib) const;

        std::string toString() const;
        bool operator==(const MutationGraph&) const;

};
