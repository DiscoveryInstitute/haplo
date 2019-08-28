#pragma once
#include "Blocks.h"

using LocusId = long;
using BlockId = long;
using AlleleId = long;
using Chromosome = Blocks<AlleleId>;
using ChromosomeId = long;

static constexpr AlleleId UNKNOWN_ALLELE{};
static constexpr long UNKNOWN_GENERATION = -1;

enum ChromosomeType { AUTO, X, Y, MITO };

