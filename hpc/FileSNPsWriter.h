#pragma once
#include <random>
#include <vector>
#include <utility>
#include "core/HaploBlockVariants.h"
#include "core/PersonNode.h"

class FileSNPsWriter {
public:
    static void WriteSNPs(
            const std::string& output_filename,
            const HaploBlockVariants& variants,
            const std::pair<std::vector<PersonNode>, std::vector<PersonNode>>& whole_generation,
            std::mt19937_64& rng);
};

