#pragma once
#include <initializer_list>
#include <map>
#include <random>
#include <set>
#include <vector>
#include "BlocksRange.h"
#include "Chromosome.h"

/**
    HaploBlockBoundaries:
    Stores the boundaries where recombination events could happen, or are recorded to have
    happened, along the simulated chromosome.
    Used to generate new random recombination patterns and primordial haplotype block patterns
    (a haplotype block is a section of chromosome within which no recombination happens).
    Records which boundaries are actually used in order to define observable haplotype blocks
    (sections of chromosomes where no ancestral recombinations were observed to happen even if
    they logically could have).
*/
class HaploBlockBoundaries {

        friend class FileCheckpointer;

    private:

        std::vector<LocusId> m_block_boundaries;
        std::map<LocusId,long> m_block_boundary_used;

	public:

        HaploBlockBoundaries() {} // for creating unassigned variables

        HaploBlockBoundaries(long chromosome_length, long max_blocks, std::mt19937_64& rng);

        BlocksRange createRandomRecombinationPattern(double rate, std::mt19937_64& rng) const;

        long getChromosomeLength() const { return m_block_boundaries.back(); }
        void recordBlockBoundaries(const BlocksRange&);
        std::vector<LocusId> getRecordedBlockBoundaries() const;
        std::vector<LocusId> getRandomBlockBoundaries(long nblocks, std::mt19937_64& rng) const;

		std::string toString() const;
        bool operator==(const HaploBlockBoundaries&) const;

};
