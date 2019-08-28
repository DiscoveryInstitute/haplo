#pragma once
#include <ctime>
#include <random>
#include <string>
#include "core/AncestralChromosomes.h"
#include "core/AncestralRecombinationGraph.h"
#include "core/HaploBlockBoundaries.h"
#include "core/HaploBlockVariants.h"
#include "core/MutationRateHistory.h"
#include "core/Parser.h"
#include "core/PopulationStructureHistory.h"

class FileCheckpointer {

    private:

        const std::string m_filename;
        const long minimum_interval;
        long previous_time;

    public:

        FileCheckpointer(const std::string& filename, long min_int)
            : m_filename(filename),
              minimum_interval(min_int*CLOCKS_PER_SEC),
              previous_time(-min_int*CLOCKS_PER_SEC) {}

        bool timeDue();

        void save( const bool& forwards, const long& step,
                   const std::mt19937_64& rng,
                   const PopulationStructureHistory& psh,
                   const MutationRateHistory& mrh,
                   const HaploBlockBoundaries& hbb,
                   const AncestralRecombinationGraph& arg,
                   const HaploBlockVariants& hbv,
                   const AncestralChromosomes& ac);

        void load( bool& forwards, long& step,
                   std::mt19937_64& rng,
                   PopulationStructureHistory& psh,
                   MutationRateHistory& mrh,
                   HaploBlockBoundaries& hbb,
                   AncestralRecombinationGraph& arg,
                   HaploBlockVariants& hbv,
                   AncestralChromosomes& ac);

};
