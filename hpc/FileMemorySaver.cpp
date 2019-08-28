#include "FileMemorySaver.h"

#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <vector>
#include "core/Blocks.h"
#include "core/BlocksRange.h"
#include "FileDataPrimitiveIO.h"

using namespace std;

constexpr unsigned char FATHER_ID_DUMMY = 1<<0;
constexpr unsigned char MOTHER_ID_DUMMY = 1<<1;
constexpr unsigned char FATHER_REC_EMPTY = 1<<2;
constexpr unsigned char FATHER_REC_FULL  = 1<<3;
constexpr unsigned char MOTHER_REC_EMPTY = 1<<4;
constexpr unsigned char MOTHER_REC_FULL  = 1<<5;
constexpr unsigned char FATHER_EXT_EMPTY = 1<<6;
constexpr unsigned char MOTHER_EXT_EMPTY = 1<<7;
constexpr long DUMMY_ID = AncestralRecombinationGraph::DUMMY_ID;
#define NOOP /* meaning: no-operation; do nothing */

FileMemorySaver::FileMemorySaver(const string& pfx) : m_file_prefix(pfx) {
    ofstream ofs(pfx + "test");
    if (!ofs) throw runtime_error("Directory not available to write files: " + pfx);
}


void FileMemorySaver::unloadGenerationToFile(AncestralRecombinationGraph& arg, long gen) const {
    ofstream ofs(m_file_prefix+"ARG"+to_string(gen)+".bin", ios::out | ios::binary);
    auto& males = arg.m_nodes[gen].first;
    auto& females = arg.m_nodes[gen].second;
    auto nmales = males.size();
    auto nfemales = females.size();
    writePrimitive(ofs, nmales);
    writePrimitive(ofs, nfemales);
    for (auto* persons : { &males, &females } )
    for (auto& person : *persons) {
        const auto& fatherId = person.fatherId;
        const auto& motherId = person.motherId;
        const auto& fatherRec = person.chromosomeFromFatherRecombinationPattern;
        const auto& motherRec = person.chromosomeFromMotherRecombinationPattern;
        const auto& fatherExt = person.chromosomeFromFatherExtancyPattern;
        const auto& motherExt = person.chromosomeFromMotherExtancyPattern;
        unsigned char flags = 0;
        if (fatherId==DUMMY_ID) flags |= FATHER_ID_DUMMY;
        if (motherId==DUMMY_ID) flags |= MOTHER_ID_DUMMY;
        if (fatherRec.empty()) flags |= FATHER_REC_EMPTY;
        if (fatherRec.full())  flags |= FATHER_REC_FULL;
        if (motherRec.empty()) flags |= MOTHER_REC_EMPTY;
        if (motherRec.full())  flags |= MOTHER_REC_FULL;
        if (fatherExt.empty()) flags |= FATHER_EXT_EMPTY;
        if (motherExt.empty()) flags |= MOTHER_EXT_EMPTY;
        writePrimitive(ofs, flags);
        if (flags & FATHER_ID_DUMMY) NOOP;
        else                         writePrimitive(ofs, fatherId);
        if (flags & MOTHER_ID_DUMMY) NOOP;
        else                         writePrimitive(ofs, motherId);
        if      (flags & FATHER_REC_EMPTY) NOOP;
        else if (flags & FATHER_REC_FULL)  NOOP;
        else                               writePrimitiveVector(ofs, fatherRec.m_toggles);
        if      (flags & MOTHER_REC_EMPTY) NOOP;
        else if (flags & MOTHER_REC_FULL)  NOOP;
        else                               writePrimitiveVector(ofs, motherRec.m_toggles);
        if (flags & FATHER_EXT_EMPTY) NOOP;
        else                          writePrimitiveVector(ofs, fatherExt.m_toggles);
        if (flags & MOTHER_EXT_EMPTY) NOOP;
        else                          writePrimitiveVector(ofs, motherExt.m_toggles);
    }
    ofs.close();
    typedef AncestralRecombinationGraph::ChildNode T;
    males = vector<T>();
    females = vector<T>();
}


void FileMemorySaver::reloadGenerationFromFile(AncestralRecombinationGraph& arg, long gen) const {
    ifstream ifs(m_file_prefix+"ARG"+to_string(gen)+".bin", ios::in | ios::binary);
    long nmales, nfemales;
    readPrimitive(ifs, nmales);
    readPrimitive(ifs, nfemales);
    auto& males = arg.m_nodes[gen].first;
    auto& females = arg.m_nodes[gen].second;
    males.resize(nmales);
    females.resize(nfemales);
    for (auto* persons : { &males, &females } )
    for (auto& person : *persons) {
        auto& fatherId = person.fatherId;
        auto& motherId = person.motherId;
        auto& fatherRec = person.chromosomeFromFatherRecombinationPattern;
        auto& motherRec = person.chromosomeFromMotherRecombinationPattern;
        auto& fatherExt = person.chromosomeFromFatherExtancyPattern;
        auto& motherExt = person.chromosomeFromMotherExtancyPattern;
        unsigned char flags;
        readPrimitive(ifs, flags);
        if (flags & FATHER_ID_DUMMY) fatherId = DUMMY_ID;
        else                         readPrimitive(ifs, person.fatherId);
        if (flags & MOTHER_ID_DUMMY) motherId = DUMMY_ID;
        else                         readPrimitive(ifs, person.motherId);
        if      (flags & FATHER_REC_EMPTY) fatherRec = BlocksRange::createEmpty();
        else if (flags & FATHER_REC_FULL)  fatherRec = BlocksRange::createFull();
        else                               readPrimitiveVector(ifs, fatherRec.m_toggles);
        if      (flags & MOTHER_REC_EMPTY) motherRec = BlocksRange::createEmpty();
        else if (flags & MOTHER_REC_FULL)  motherRec = BlocksRange::createFull();
        else                               readPrimitiveVector(ifs, motherRec.m_toggles);
        if (flags & FATHER_EXT_EMPTY) fatherExt = BlocksRange::createEmpty();
        else                          readPrimitiveVector(ifs, fatherExt.m_toggles);
        if (flags & MOTHER_EXT_EMPTY) motherExt = BlocksRange::createEmpty();
        else                          readPrimitiveVector(ifs, motherExt.m_toggles);
    }
}


void FileMemorySaver::discardGeneration(AncestralRecombinationGraph& arg, long gen) {
    typedef AncestralRecombinationGraph::ChildNode T;
    arg.m_nodes[gen] = make_pair(vector<T>(), vector<T>());
    m_discard_list.push_back(m_file_prefix+"ARG"+to_string(gen)+".bin");
}

void FileMemorySaver::unloadGenerationToFile(AncestralChromosomes& ac, long gen) const {
    ofstream ofs(m_file_prefix+"AC"+to_string(gen)+".bin", ios::out | ios::binary);
    auto& males = ac.m_nodes[gen].first;
    auto& females = ac.m_nodes[gen].second;
    auto nmales = males.size();
    auto nfemales = females.size();
    writePrimitive(ofs, nmales);
    writePrimitive(ofs, nfemales);
    for (auto* persons : { &males, &females } )
    for (auto& person : *persons) {
        writePrimitiveVector(ofs, person.chromosomeFromFather.m_intervals);
        writePrimitiveVector(ofs, person.chromosomeFromMother.m_intervals);
    }
    ofs.close();
    typedef PersonNode T;
    males = vector<T>();
    females = vector<T>();
}

void FileMemorySaver::reloadGenerationFromFile(AncestralChromosomes& ac, long gen) const {
    ifstream ifs(m_file_prefix+"AC"+to_string(gen)+".bin", ios::in | ios::binary);
    long nmales, nfemales;
    readPrimitive(ifs, nmales);
    readPrimitive(ifs, nfemales);
    auto& males = ac.m_nodes[gen].first;
    auto& females = ac.m_nodes[gen].second;
    males.resize(nmales);
    females.resize(nfemales);
    for (auto* persons : { &males, &females } )
    for (auto& person : *persons) {
        readPrimitiveVector(ifs, person.chromosomeFromFather.m_intervals);
        readPrimitiveVector(ifs, person.chromosomeFromMother.m_intervals);
    }
}

void FileMemorySaver::discardGeneration(AncestralChromosomes& ac, long gen) {
    typedef PersonNode T;
    ac.m_nodes[gen] = make_pair(vector<T>(), vector<T>());
    m_discard_list.push_back(m_file_prefix+"AC"+to_string(gen)+".bin");
}

void FileMemorySaver::deleteDiscardedGenerations() {
    for (string& filename : m_discard_list) remove(filename.c_str());
    m_discard_list.clear();
}

