#pragma once
#include <utility>
#include <vector>
#include "GraphDataBinned.h"
#include "GraphDataBinned2D.h"
#include "HaploBlockVariants.h"
#include "Logger.h"
#include "Parser.h"
#include "PersonNode.h"

class Analysis {

    public:

    long total_blocks;                          // Q in the paper
    long total_snp_sites;                       // S in the paper
    double sn_diversity;                        // Pi in the paper
    GraphDataBinned snp_freq_spectrum;                    // Phi in the paper
    GraphDataBinned snp_distribution;                     // Not in paper
    GraphDataBinned correlation_function;                 // r^2 in the paper
    GraphDataBinned dprime_function;                      // D' in the paper
    GraphDataBinned complete_linkage_function;            // CLDP in the paper
    GraphDataBinned almost_complete_linkage_function;     // Not in paper
    GraphDataBinned normalized_kullback_leibler_function; // Not in paper
    GraphDataBinned sigma_squared_function;               // Not in paper
    GraphDataBinned2D freq_pair_spectrum;        // Not in paper

    static Analysis Calculate(
            Logger& logger,
            Parser& inputs,
            HaploBlockVariants& variants,
            std::pair<std::vector<PersonNode>, std::vector<PersonNode>>& whole_gen);

};

