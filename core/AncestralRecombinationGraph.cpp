#include "AncestralRecombinationGraph.h"
#include <memory>
#include <random>
#include <utility>
#include <vector>
#include "BlocksRange.h"
#include "ParentPairUrn.h"
#include "PopulationStructureHistory.h"

using namespace std;


template<typename T> vector<T> trim(vector<T>& v, long n) {
	return vector<T>(make_move_iterator(v.begin()), make_move_iterator(v.begin() + n)); }

static constexpr long DID = AncestralRecombinationGraph::DUMMY_ID;
bool AncestralRecombinationGraph::ChildNode::extantFather() const {
    return !chromosomeFromFatherExtancyPattern.empty();
}
bool AncestralRecombinationGraph::ChildNode::extantMother() const {
    return !chromosomeFromMotherExtancyPattern.empty();
}
bool AncestralRecombinationGraph::ChildNode::extant() const {
    return extantFather() || extantMother();
}
bool AncestralRecombinationGraph::ChildNode::operator==(const ChildNode& o) const {
    return fatherId == o.fatherId
        && motherId == o.motherId
        && chromosomeFromFatherRecombinationPattern == o.chromosomeFromFatherRecombinationPattern
		&& chromosomeFromMotherRecombinationPattern == o.chromosomeFromMotherRecombinationPattern
		&& chromosomeFromFatherExtancyPattern == o.chromosomeFromFatherExtancyPattern
		&& chromosomeFromMotherExtancyPattern == o.chromosomeFromMotherExtancyPattern;
}


/**
 *  Initialize children of the extant (final) generation.
 */
AncestralRecombinationGraph& AncestralRecombinationGraph::initialiseExtantGeneration(
        const pair<long,long> samples)
{

    // First check the requested sample is consistent with population size.
    const long extantGeneration = 0L;
    long nMalesSample = samples.first;
    long nFemalesSample = samples.second;
    long nMalesPopn = m_history->getNumberOfMales(extantGeneration);
    long nFemalesPopn = m_history->getNumberOfFemales(extantGeneration);
    if (nMalesSample > nMalesPopn || nFemalesSample > nFemalesPopn) {
        throw invalid_argument("Sample must be no larger than extant population.");
    }

    // Specify which chromosomes are of interest or to be ignored for each type.
    //   Relevant chromosomes are by definition wholly extant (ALL).
    //   Ignored chromosomes are treated as if wholly non-extant (NONE).
    const BlocksRange ALL = BlocksRange::createFull();
    const BlocksRange NONE = BlocksRange();
    const BlocksRange& DRP = NONE; // DUMMY_RECOMBINATION_PATTERN
    ChildNode MALE_NODE;
    ChildNode FEMALE_NODE;
    switch (m_chromosome_type) {
    case AUTO: // Both male and female have autosome from both father and mother.
        MALE_NODE   = { DID, DID, DRP, DRP, ALL, ALL };
        FEMALE_NODE = { DID, DID, DRP, DRP, ALL, ALL };
        break;
    case X:    // Male has X chromosome from mother. Female has X from both father and mother.
        MALE_NODE   = { DID, DID, DRP, DRP, NONE, ALL };
        FEMALE_NODE = { DID, DID, DRP, DRP, ALL,  ALL };
        break;
    case Y:    // Male has Y chromosome from father. Female has no Y chromosome.
        MALE_NODE   = { DID, DID, DRP, DRP, ALL,  NONE };
        FEMALE_NODE = { DID, DID, DRP, DRP, NONE, NONE };
        break;
    case MITO: // Both male and female have DNA from the mother only.
        MALE_NODE   = { DID, DID, DRP, DRP, NONE, ALL };
        FEMALE_NODE = { DID, DID, DRP, DRP, NONE, ALL };
        break;
    }

    // Create the representation of the final/extant/sampled generation.
    auto maleChildNodes = vector<ChildNode>(nMalesSample, MALE_NODE);
    auto femaleChildNodes = vector<ChildNode>(nFemalesSample, FEMALE_NODE);
    m_nodes[extantGeneration] = { move(maleChildNodes), move(femaleChildNodes) };
    return *this;
}

AncestralRecombinationGraph& AncestralRecombinationGraph::backwardsCoalesceAllGenerations(
        const double alpha, const double beta, const double rate,
        mt19937_64& rng,
        const bool cull_non_parents,
        const bool cull_non_ancestral,
        const bool hide_non_ancestral) {

    const long nGenerations = m_history->getNumberOfGenerations();
    const long extantGeneration = 0L;
    const long foundingGeneration = nGenerations-1;
    for (long gen=extantGeneration; gen<foundingGeneration; ++gen) {
        backwardsCoalesceOneGeneration(
                    gen, alpha, beta, rate,
                    rng, cull_non_parents, cull_non_ancestral, hide_non_ancestral);
    }
    return *this;
}

AncestralRecombinationGraph& AncestralRecombinationGraph::backwardsCoalesceOneGeneration(
        const long gen,
        const double alpha, const double beta, const double rate,
        mt19937_64& rng,
        const bool cull_non_parents,
        const bool cull_non_ancestral,
        const bool hide_non_ancestral) {

    const BlocksRange ALL = BlocksRange::createFull();
    const BlocksRange NONE = BlocksRange();
    const BlocksRange& INHERIT_GRANDFATHER = ALL;  // (all of father.chromosomeFromFather)
    const BlocksRange& INHERIT_GRANDMOTHER = NONE; // (none of father.chromosomeFromFather)
    const BlocksRange& DRP = NONE; // DUMMY_RECOMBINATION_PATTERN
    const ChildNode NON_EXTANT_NODE { DID, DID, DRP, DRP, NONE, NONE };

    const long nGenerations = m_history->getNumberOfGenerations();
    const long lastGeneration = nGenerations-1;
    HaploBlockBoundaries* hbb = m_haplo_block_boundaries;
    if (gen < 0L) throw invalid_argument("Gen must be >=0");
    if (gen >= lastGeneration) throw invalid_argument("Gen must have children in ARG");


    auto& maleChildNodes = m_nodes[gen].first;
    auto& femaleChildNodes = m_nodes[gen].second;

    // Adults are potential parents. Number of actual parents is generated in the parent urn.
    // All parent chromosomes are completely non-extant by default,
    // (unless they have children inherited extant blocks from them - see below).
    long nMaleAdults = m_history->getNumberOfMales(gen+1);
    long nFemaleAdults = m_history->getNumberOfFemales(gen+1);
    auto fatherNodes = vector<ChildNode>(nMaleAdults, NON_EXTANT_NODE);
    auto motherNodes = vector<ChildNode>(nFemaleAdults, NON_EXTANT_NODE);

    ParentPairUrn parentUrn(nMaleAdults, nFemaleAdults, alpha, beta);

    for (auto* childNodes : { &maleChildNodes, &femaleChildNodes })
    for (auto& childNode : *childNodes) {

        auto parentIds = parentUrn.pick(rng);
        long fatherId = parentIds.first;
        long motherId = parentIds.second;
        auto& fatherEP = childNode.chromosomeFromFatherExtancyPattern;
        auto& motherEP = childNode.chromosomeFromMotherExtancyPattern;
        bool fatherCPE = !fatherEP.empty(); // father chromosome partially extant
        bool motherCPE = !motherEP.empty(); // mother chromosome partially extant
        BlocksRange fatherRecombinations;
        BlocksRange motherRecombinations;
        switch (m_chromosome_type) {
        case AUTO: // Recombination on chromosome from father and from mother
            if (fatherCPE) fatherRecombinations = hbb->createRandomRecombinationPattern(rate, rng);
            if (motherCPE) motherRecombinations = hbb->createRandomRecombinationPattern(rate, rng);
            break;
        case X:    // Recombination from mother; from father is always from grandmother.
            if (fatherCPE) fatherRecombinations = INHERIT_GRANDMOTHER;
            if (motherCPE) motherRecombinations = hbb->createRandomRecombinationPattern(rate, rng);
            break;
        case Y:    // From father is always from grandfather. No inheritance from mother.
            if (fatherCPE) fatherRecombinations = INHERIT_GRANDFATHER;
            if (motherCPE) motherRecombinations = DRP;
            break;
        case MITO: // From mother is always from grandmother. No inheritance from father.
            if (fatherCPE) fatherRecombinations = DRP;
            if (motherCPE) motherRecombinations = INHERIT_GRANDMOTHER;
            break;
        }

        if (fatherCPE) {
            auto& fathersFatherEP = fatherNodes[fatherId].chromosomeFromFatherExtancyPattern;
            auto& fathersMotherEP = fatherNodes[fatherId].chromosomeFromMotherExtancyPattern;
            auto fathersFatherPartialEP = fatherEP.intersectionWith(fatherRecombinations);
            auto fathersMotherPartialEP = fatherEP.intersectionWithInverse(fatherRecombinations);
            if (hide_non_ancestral) {
                if (fathersFatherPartialEP.empty()) fatherRecombinations = INHERIT_GRANDMOTHER;
                if (fathersMotherPartialEP.empty()) fatherRecombinations = INHERIT_GRANDFATHER;
            }
            fathersFatherEP = fathersFatherEP.unionWith(fathersFatherPartialEP);
            fathersMotherEP = fathersMotherEP.unionWith(fathersMotherPartialEP);
            hbb->recordBlockBoundaries(fatherRecombinations);
        }
        if (motherCPE) {
            auto& mothersFatherEP = motherNodes[motherId].chromosomeFromFatherExtancyPattern;
            auto& mothersMotherEP = motherNodes[motherId].chromosomeFromMotherExtancyPattern;
            auto mothersFatherPartialEP = motherEP.intersectionWith(motherRecombinations);
            auto mothersMotherPartialEP = motherEP.intersectionWithInverse(motherRecombinations);
            if (hide_non_ancestral) {
                if (mothersFatherPartialEP.empty()) motherRecombinations = INHERIT_GRANDMOTHER;
                if (mothersMotherPartialEP.empty()) motherRecombinations = INHERIT_GRANDFATHER;
            }
            mothersFatherEP = mothersFatherEP.unionWith(mothersFatherPartialEP);
            mothersMotherEP = mothersMotherEP.unionWith(mothersMotherPartialEP);
            hbb->recordBlockBoundaries(motherRecombinations);
        }

        childNode.fatherId = fatherId;
        childNode.motherId = motherId;
        childNode.chromosomeFromFatherRecombinationPattern = move(fatherRecombinations);
        childNode.chromosomeFromMotherRecombinationPattern = move(motherRecombinations);

    }

    // Clean-up
    long nFathers = parentUrn.getNumberOfFathers();
    long nMothers = parentUrn.getNumberOfMothers();
    if (cull_non_ancestral) {
        // Shift all the ancestral parents down the vector/array,
        // overwriting any non-ancestral parents; keep track of id changes.
        vector<long> fatherIdNew(nFathers);
        vector<long> motherIdNew(nMothers);
        long fidNew = 0, midNew = 0;
        for (long id=0; id<nFathers; id++) {
            if (fatherNodes[id].extant()) {
                fatherNodes[fidNew] = fatherNodes[id];
                fatherIdNew[id] = fidNew++;
            } else {
                fatherIdNew[id] = DID;
            }
        }
        for (long id=0; id<nMothers; id++) {
            if (motherNodes[id].extant()) {
                motherNodes[midNew] = motherNodes[id];
                motherIdNew[id] = midNew++;
            } else {
                motherIdNew[id] = DID;
            }
        }
        // Update the parent id changes in the children.
        for (auto* childNodes : { &maleChildNodes, &femaleChildNodes })
        for (auto& childNode : *childNodes) {
            childNode.fatherId = fatherIdNew[childNode.fatherId];
            childNode.motherId = motherIdNew[childNode.motherId];
        }
        // Trim the vector/array.
        nFathers = fidNew;
        nMothers = midNew;
        fatherNodes = trim(fatherNodes, nFathers);
        motherNodes = trim(motherNodes, nMothers);
    }
    else if (cull_non_parents) {
        fatherNodes = trim(fatherNodes, nFathers);
        motherNodes = trim(motherNodes, nMothers);
    }
    if (!hide_non_ancestral) {
        for (auto* nodes : { &fatherNodes, &motherNodes })
        for (auto& node : *nodes) {
            node.chromosomeFromFatherExtancyPattern = ALL;
            node.chromosomeFromMotherExtancyPattern = ALL;
        }
    }

    m_nodes[gen+1] = { fatherNodes, motherNodes };

    return *this;
}



long AncestralRecombinationGraph::getNumberOfMalesWithAncestralHaploBlocks(long gen) const {
    auto& males = m_nodes[gen].first;
    long sum = 0L;
    for(auto& male : males) {
        if (male.extant()) ++sum;
    }
    return sum;
}
long AncestralRecombinationGraph::getNumberOfFemalesWithAncestralHaploBlocks(long gen) const {
    auto& females = m_nodes[gen].second;
    long sum = 0L;
    for(auto& female : females) {
        if (female.extant()) ++sum;
    }
    return sum;
}
long AncestralRecombinationGraph::getSizeOfAncestralHaploBlocks(long gen) const {
    const long max_range = m_haplo_block_boundaries->getChromosomeLength();
    long sum = 0L;
    auto& males = m_nodes[gen].first;
    for(auto& male : males) {
        sum += male.chromosomeFromFatherExtancyPattern.getSizeOfRange(max_range);
        sum += male.chromosomeFromMotherExtancyPattern.getSizeOfRange(max_range);
    }
    auto& females = m_nodes[gen].second;
    for(auto& female : females) {
        sum += female.chromosomeFromFatherExtancyPattern.getSizeOfRange(max_range);
        sum += female.chromosomeFromMotherExtancyPattern.getSizeOfRange(max_range);
    }
    return sum;
}

long AncestralRecombinationGraph::getTotalNumberOfMales() const {
    const long ngens = getNumberOfGenerations();
    long sum = 0L;
    for (long gen=0; gen<ngens; ++gen) sum += getMales(gen).size();
    return sum;
}
long AncestralRecombinationGraph::getTotalNumberOfFemales() const {
    const long ngens = getNumberOfGenerations();
    long sum = 0L;
    for (long gen=0; gen<ngens; ++gen) sum += getFemales(gen).size();
    return sum;
}
long AncestralRecombinationGraph::getTotalNumberOfMalesWithAncestralHaploBlocks() const {
    const long ngens = getNumberOfGenerations();
    long sum = 0L;
    for (long gen=0; gen<ngens; ++gen) sum += getNumberOfMalesWithAncestralHaploBlocks(gen);
    return sum;
}
long AncestralRecombinationGraph::getTotalNumberOfFemalesWithAncestralHaploBlocks() const {
    const long ngens = getNumberOfGenerations();
    long sum = 0L;
    for (long gen=0; gen<ngens; ++gen) sum += getNumberOfFemalesWithAncestralHaploBlocks(gen);
    return sum;
}
long AncestralRecombinationGraph::getTotalSizeOfAncestralHaploBlocks() const {
    const long ngens = getNumberOfGenerations();
    long sum = 0L;
    for (long gen=0; gen<ngens; ++gen) sum += getSizeOfAncestralHaploBlocks(gen);
    return sum;
}

size_t AncestralRecombinationGraph::getSizeInMemory(long gen) const {
    size_t sum = sizeof(m_nodes[gen]);
    for (auto* persons : { &m_nodes[gen].first, &m_nodes[gen].second })
    for (auto& person : *persons) {
        sum += sizeof(person.fatherId);
        sum += sizeof(person.motherId);
        sum += person.chromosomeFromFatherRecombinationPattern.getSizeInMemory();
        sum += person.chromosomeFromMotherRecombinationPattern.getSizeInMemory();
        sum += person.chromosomeFromFatherExtancyPattern.getSizeInMemory();
        sum += person.chromosomeFromMotherExtancyPattern.getSizeInMemory();
    }
    return sum;
}

bool AncestralRecombinationGraph::operator==(const AncestralRecombinationGraph& o) const {
    return ( m_history == o.m_history
             || ( m_history != nullptr
                && o.m_history != nullptr
                && *m_history == *o.m_history) )
        && ( m_haplo_block_boundaries == o.m_haplo_block_boundaries
             || ( m_haplo_block_boundaries != nullptr
                && o.m_haplo_block_boundaries != nullptr
                && *m_haplo_block_boundaries == *o.m_haplo_block_boundaries) )
        && m_nodes           ==  o.m_nodes
        && m_chromosome_type ==  o.m_chromosome_type;
}
