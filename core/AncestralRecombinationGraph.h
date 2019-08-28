#pragma once

#include <cstddef>
#include <utility>
#include <vector>
#include "BlocksRange.h"
#include "HaploBlockBoundaries.h"
#include "PopulationStructureHistory.h"

/**
    AncestralRecombinationGraph:
    Stores the Ancestral Recombination Graph for the whole population by storing, for each
    ancestral individual, the ids of their parents in the previous generation,
    the recombination pattern by which the chromosome inherited from the father was recombined
    from the father's two chromosomes, and a mask showing which parts of that chromosome are
    *ancestral*: that is, which haplotype blocks have direct direct descendants in any
    chromosomes of any person in the final sample, and the corresponding recombination pattern
    and mask for the chromosome inherited from the mother.
    The ARG is built backwards in time from the final generation (this is possible because
    the recombination process is symmetric in time, unlike a mutation or selection process).
    It can also be used to query how much genetic information present at a given timepoint is
    ancestral to the final sample.
*/
class AncestralRecombinationGraph {

        friend class FileCheckpointer;
        friend class FileMemorySaver;

    public:

		/**
		 *  ChildNode contains the minimum information required to generate the chromosomes of a
		 *  child from the previous generation.
		 *
		 *  The final AncestralRecombinationGraph is a directed network of ChildNodes.
		 */
		struct ChildNode {
			long fatherId;
			long motherId;
			BlocksRange chromosomeFromFatherRecombinationPattern;
			BlocksRange chromosomeFromMotherRecombinationPattern;
			BlocksRange chromosomeFromFatherExtancyPattern;
			BlocksRange chromosomeFromMotherExtancyPattern;
            bool extantFather() const;
            bool extantMother() const;
            bool extant() const;
            bool operator==(const ChildNode&) const;
        };

	private:

		/** The sizes and structures of populations through time. */
        PopulationStructureHistory* m_history;

        /** The actual haplotype boundaries observed (randomly created). */
        HaploBlockBoundaries* m_haplo_block_boundaries;

        /** The nodes of the graph. */
        std::vector<std::pair<std::vector<ChildNode>, std::vector<ChildNode>>> m_nodes;

        /** Autosome, X, Y or mitochondrial.  */
        ChromosomeType m_chromosome_type;

  public:

        AncestralRecombinationGraph() // for creating unassigned variables
            : m_history(nullptr), m_haplo_block_boundaries(nullptr), m_chromosome_type(AUTO) {}

        /**
         *  Constructor for AncestralRecombinationGraph with given inputs:
         *  population size history, and maximum number of haplotype blocks.
         */
        AncestralRecombinationGraph(
                PopulationStructureHistory* history,
                HaploBlockBoundaries* boundaries,
                ChromosomeType chromosome_type = AUTO)
            : m_history(history),
              m_haplo_block_boundaries(boundaries),
              m_nodes(history->getNumberOfGenerations()),
              m_chromosome_type(chromosome_type)
            {}

        AncestralRecombinationGraph& initialiseExtantGeneration(const std::pair<long,long> samples);
        AncestralRecombinationGraph& initialiseExtantGeneration(const long samples) {
            return initialiseExtantGeneration({samples/2, samples/2});
        }
        /**
         *  Randomly generate the AncestralRecombinationGraph given these additional
         *  input parameters which can be specified per generation:
		 *  'alpha' influences the number of females that become mothers.
		 *  'beta' influences the number of males that become fathers.
         *  'recombination_rate' is the expected number of recombinations per copied chromosome.
		 */
        AncestralRecombinationGraph& backwardsCoalesceAllGenerations(
                const double alpha, const double beta,
                const double recombination_rate,
                std::mt19937_64& rng,
                const bool cull_non_parents = true,
                const bool cull_non_ancesral = true,
                const bool hide_non_ancestral = true);
        AncestralRecombinationGraph& backwardsCoalesceOneGeneration(
                const long gen,
                const double alpha, const double beta,
                const double recombination_rate,
                std::mt19937_64& rng,
                const bool cull_non_parents = true,
                const bool cull_non_ancesral = true,
                const bool hide_non_ancestral = true);

		long getNumberOfGenerations() const { return m_nodes.size(); }
        const std::vector<ChildNode>& getMales(long gen) const { return m_nodes[gen].first; }
        const std::vector<ChildNode>& getFemales(long gen) const { return m_nodes[gen].second; }

        // Statistics methods
        long getNumberOfMalesWithAncestralHaploBlocks(long gen) const;
        long getNumberOfFemalesWithAncestralHaploBlocks(long gen) const;
        long getSizeOfAncestralHaploBlocks(long gen) const;

        long getTotalNumberOfMales() const;
        long getTotalNumberOfFemales() const;
        long getTotalNumberOfMalesWithAncestralHaploBlocks() const;
        long getTotalNumberOfFemalesWithAncestralHaploBlocks() const;
        long getTotalSizeOfAncestralHaploBlocks() const;

        size_t getSizeInMemory(long gen) const;
        bool operator==(const AncestralRecombinationGraph&) const;

        static constexpr long DUMMY_ID = -1L;

};
