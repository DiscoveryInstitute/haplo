#pragma once
#include "Chromosome.h"

struct PersonNode {
    Chromosome chromosomeFromFather;
    Chromosome chromosomeFromMother;

    PersonNode() {}
    PersonNode(Chromosome&& f, Chromosome&& m)
        : chromosomeFromFather(f), chromosomeFromMother(m) {}

    bool operator==(const PersonNode& o) const {
        return chromosomeFromFather == o.chromosomeFromFather
            && chromosomeFromMother == o.chromosomeFromMother;
    }
};
