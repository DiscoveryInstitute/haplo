#pragma once

#include <cstddef>
#include "AncestralRecombinationGraph.h"
#include "Chromosome.h"
#include "HaploBlockVariants.h"
#include "PersonNode.h"
#include "PopulationStructureHistory.h"

/**
    AncestralChromosomes:
    Generates and stores the ancestral genetic material of the whole population through all
    generations by storing, for each ancestral individual, a representation of their two
    chromosomes, but masked by the 'extancy patterns' that were calculated in the
    AncestralRecombinationGraph, in order to reduce the amount of information stored.
    Chromosomes are represented by a set of sequential pairs: block-locus and the allele at that
    locus. Masked / non-extant regions of chromosomes are represented by UNKNOWN_ALLELE.
    Contiguous regions with the same allele id are merged into one block-loci/allele-id pair.
    In practice, the large contiguous regions of UNKNOWN_ALLELE allow the representation to be
    greatly compressed in memory.
*/
class AncestralChromosomes {

        friend class FileCheckpointer;
        friend class FileMemorySaver;

    private:

        AncestralRecombinationGraph* m_graph;
        HaploBlockVariants* m_variants;

        // The nodes of the graph.
		std::vector<std::pair<std::vector<PersonNode>, std::vector<PersonNode>>> m_nodes;

	public:

        AncestralChromosomes() : m_graph(nullptr), m_variants(nullptr) {} // for creating unassigned variables
        bool unassigned() const { return m_nodes.empty(); }

        AncestralChromosomes(
                AncestralRecombinationGraph* graph,
                HaploBlockVariants* variants)
            : m_graph(graph), m_variants(variants), m_nodes(graph->getNumberOfGenerations())
            {}

        AncestralChromosomes& initialiseFoundingGeneration();
        AncestralChromosomes& forwardsDropAllGenerations(
                const double mu, std::mt19937_64& rng);
        AncestralChromosomes& forwardsDropOneGeneration(const long gen,
                const double mu, std::mt19937_64& rng);

        long getNumberOfGenerations() const { return m_nodes.size(); }
        const std::vector<PersonNode>& getMales(long gen) const { return m_nodes[gen].first; }
        const std::vector<PersonNode>& getFemales(long gen) const { return m_nodes[gen].second; }

        const std::pair<std::vector<PersonNode>, std::vector<PersonNode>>&
            getGeneration(long gen) const { return m_nodes[gen]; }

        std::pair<std::vector<PersonNode>, std::vector<PersonNode>>&
            getExtantGeneration() { return m_nodes[0]; }

        size_t getSizeInMemory(long gen) const;

        std::string toString(long gen) const;
        bool operator==(const AncestralChromosomes&) const;

};
