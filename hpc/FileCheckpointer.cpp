#include "FileCheckpointer.h"
#include <algorithm>
#include <type_traits>
#define typeof(x) std::remove_reference<decltype((x))>::type
#include "FileDataPrimitiveIO.h"

using namespace std;

bool FileCheckpointer::timeDue() {
    long current_time = clock();
    return (current_time >= previous_time + minimum_interval);
}

void FileCheckpointer::save(const bool& forwards, const long& step,
                            const mt19937_64& rng,
                            const PopulationStructureHistory& psh,
                            const MutationRateHistory& mrh,
                            const HaploBlockBoundaries& hbb,
                            const AncestralRecombinationGraph& arg,
                            const HaploBlockVariants& hbv,
                            const AncestralChromosomes& ac) {

    long current_time = clock();
    previous_time = current_time;

    ofstream ofs(m_filename, ios::out | ios::binary);
    writePrimitive(ofs, forwards);
    writePrimitive(ofs, step);
    ofs << rng;
    /** Parser
    writePrimitive(ofs, parser.m_map.size());
    for (auto& entry : parser.m_map) {
        writeString(ofs, entry.first);
        writeString(ofs, entry.second);
    }
    writePrimitive(ofs, parser.m_used.size());
    for (auto& entry : parser.m_used) {
        writeString(ofs, entry);
    }
    */
    /** PopulationStructureHistory */
    writePrimitiveVector(ofs, psh.m_populations);
    /** MutationRateHistory */
    writePrimitiveVector(ofs, mrh.m_rates);
    /** HaploBlockBoundaries */
    writePrimitiveVector(ofs, hbb.m_block_boundaries);
    writePrimitiveMap(ofs, hbb.m_block_boundary_used);
    /** Ancestral Recombination Graph */
    {
        writePrimitive(ofs, arg.m_chromosome_type);
        const auto& gens = arg.m_nodes;
        const decltype(arg.m_nodes)::value_type nullgen;
        size_t ngens = gens.size();
        writePrimitive(ofs, ngens);
        //size_t n2write = count_if(gens.begin(), gens.end(), [](const typeof(gens)::value_type& g){ return !g.first.empty() || !g.second.empty();});
        size_t n2write = ngens - count(gens.begin(), gens.end(), nullgen);
        writePrimitive(ofs, n2write);
        for (size_t igen=0; igen<ngens; ++igen) {
            const auto& gen = gens[igen];
            if (n2write<ngens/2) {
                if(gen == nullgen) continue;
                writePrimitive(ofs, igen);
            }
            auto& males = gen.first;
            auto& females = gen.second;
            auto nmales = males.size();
            auto nfemales = females.size();
            writePrimitive(ofs, nmales);
            writePrimitive(ofs, nfemales);
            for (auto* persons : { &males, &females } )
            for (auto& person : *persons) {
                writePrimitive(ofs, person.fatherId);
                writePrimitive(ofs, person.motherId);
                writePrimitiveVector(ofs, person.chromosomeFromFatherRecombinationPattern.m_toggles);
                writePrimitiveVector(ofs, person.chromosomeFromMotherRecombinationPattern.m_toggles);
                writePrimitiveVector(ofs, person.chromosomeFromFatherExtancyPattern.m_toggles);
                writePrimitiveVector(ofs, person.chromosomeFromMotherExtancyPattern.m_toggles);
            }
        }
    }
    /** HaploBlockVariants */
    const bool have_hbv = !hbv.unassigned();
    writePrimitive(ofs, have_hbv);
    writePrimitive(ofs, hbv.m_max_nodes_per_block);
    if (have_hbv) {
        writePrimitive(ofs, hbv.m_reference_chromosomes.size());
        for (const auto& refc : hbv.m_reference_chromosomes) {
            writePrimitiveVector(ofs, refc.m_intervals);
        }
        writePrimitive(ofs, hbv.m_variants.size());
        for (const auto& mutg : hbv.m_variants) {
            writePrimitive(ofs, mutg.m_use_mutation_loci);
            writePrimitive(ofs, mutg.m_next_id);
            writePrimitive(ofs, mutg.m_nodes.size());
            for (const auto& entry : mutg.m_nodes) {
                const auto& id = entry.first;
                const auto& node = entry.second;
                writePrimitive(ofs, id);
                writePrimitive(ofs, node.parent_id);
                writePrimitive(ofs, node.n_mutations);
                writePrimitiveVector(ofs, node.mutations);
		writePrimitive(ofs, node.n_child_nodes);
                writePrimitive(ofs, node.frequency);
            }
        }
        writePrimitiveVector(ofs, hbv.m_boundaries);
        writePrimitiveMap(ofs, hbv.m_block_ids);
    }
    /** Ancestral Chromosomes */
    const bool have_ac = !ac.unassigned();
    writePrimitive(ofs, have_ac);
    if (have_ac) {
        const auto& gens = ac.m_nodes;
        const decltype(ac.m_nodes)::value_type nullgen;
        size_t ngens = gens.size();
        writePrimitive(ofs, ngens);
        size_t n2write = ngens - count(gens.begin(), gens.end(), nullgen);
        writePrimitive(ofs, n2write);
        for (size_t igen=0; igen<ngens; ++igen) {
            const auto& gen = gens[igen];
            if (n2write<ngens/2) {
                if(gen == nullgen) continue;
                writePrimitive(ofs, igen);
            }
            auto& males = gen.first;
            auto& females = gen.second;
            auto nmales = males.size();
            auto nfemales = females.size();
            writePrimitive(ofs, nmales);
            writePrimitive(ofs, nfemales);
            for (auto* persons : { &males, &females } )
            for (auto& person : *persons) {
                writePrimitiveVector(ofs, person.chromosomeFromFather.m_intervals);
                writePrimitiveVector(ofs, person.chromosomeFromMother.m_intervals);
            }
        }
    }
}


void FileCheckpointer::load(bool& forwards, long& step,
                            mt19937_64& rng,
                            PopulationStructureHistory& psh,
                            MutationRateHistory& mrh,
                            HaploBlockBoundaries& hbb,
                            AncestralRecombinationGraph& arg,
                            HaploBlockVariants& hbv,
                            AncestralChromosomes& ac) {

    ifstream ifs(m_filename, ios::in | ios::binary);
    readPrimitive(ifs, forwards);
    readPrimitive(ifs, step);
    ifs >> rng;
    /** Parser
    size_t ninputs;
    readPrimitive(ifs, ninputs);
    for (size_t i=0; i<ninputs; ++i) {
        string key;
        string value;
        readString(ifs, key);
        readString(ifs, value);
        parser.m_map[key] = value;
    }
    size_t nused;
    readPrimitive(ifs, nused);
    for (size_t i=0; i<nused; ++i) {
        string key;
        readString(ifs, key);
        parser.m_used.insert(key);
    }
    */
    /** PopulationStructureHistory */
    readPrimitiveVector(ifs, psh.m_populations);
    /** MutationRateHistory */
    readPrimitiveVector(ifs, mrh.m_rates);
    /** HaploBlockBoundaries */
    readPrimitiveVector(ifs, hbb.m_block_boundaries);
    readPrimitiveMap(ifs, hbb.m_block_boundary_used);
    /** Ancestral Recombination Graph */
    {
        arg.m_history = &psh;
        arg.m_haplo_block_boundaries = &hbb;
        readPrimitive(ifs, arg.m_chromosome_type);
        auto& gens = arg.m_nodes;
        size_t ngens;
        readPrimitive(ifs, ngens);
        gens.resize(ngens);
        size_t n2read;
        readPrimitive(ifs, n2read);
        for (size_t i=0; i<n2read; ++i) {
            size_t igen;
            if (n2read<ngens/2) readPrimitive(ifs, igen); else igen = i;
            auto& gen = gens.at(igen);
            auto& males = gen.first;
            auto& females = gen.second;
            size_t nmales;
            size_t nfemales;
            readPrimitive(ifs, nmales);
            readPrimitive(ifs, nfemales);
            males.resize(nmales);
            females.resize(nfemales);
            for (auto* persons : { &males, &females } )
            for (auto& person : *persons) {
                readPrimitive(ifs, person.fatherId);
                readPrimitive(ifs, person.motherId);
                readPrimitiveVector(ifs, person.chromosomeFromFatherRecombinationPattern.m_toggles);
                readPrimitiveVector(ifs, person.chromosomeFromMotherRecombinationPattern.m_toggles);
                readPrimitiveVector(ifs, person.chromosomeFromFatherExtancyPattern.m_toggles);
                readPrimitiveVector(ifs, person.chromosomeFromMotherExtancyPattern.m_toggles);
            }
        }
    }
    /** HaploBlockVariants */
    bool have_hbv;
    readPrimitive(ifs, have_hbv);
    readPrimitive(ifs, hbv.m_max_nodes_per_block);
    if (have_hbv) {
        size_t nrefc;
        readPrimitive(ifs, nrefc);
        hbv.m_reference_chromosomes.resize(nrefc);
        for (auto& refc : hbv.m_reference_chromosomes) {
            readPrimitiveVector(ifs, refc.m_intervals);
        }
        size_t nvar;
        readPrimitive(ifs, nvar);
        hbv.m_variants.resize(nvar);
        for (auto& mutg : hbv.m_variants) {
            readPrimitive(ifs, mutg.m_use_mutation_loci);
            readPrimitive(ifs, mutg.m_next_id);
            size_t nnodes;
            readPrimitive(ifs, nnodes);
            for (size_t i=0; i<nnodes; ++i) {
                typeof(mutg.m_nodes)::key_type id;
                readPrimitive(ifs, id);
                auto& node = mutg.m_nodes[id];
                readPrimitive(ifs, node.parent_id);
                readPrimitive(ifs, node.n_mutations);
                readPrimitiveVector(ifs, node.mutations);
                readPrimitive(ifs, node.n_child_nodes);
                readPrimitive(ifs, node.frequency);
            }
        }
        readPrimitiveVector(ifs, hbv.m_boundaries);
        readPrimitiveMap(ifs, hbv.m_block_ids);
    }
    /** Ancestral Chromosomes */
    bool have_ac;
    readPrimitive(ifs, have_ac);
    if (have_ac) {
        ac.m_graph = &arg;
        ac.m_variants = &hbv;
        auto& gens = ac.m_nodes;
        size_t ngens;
        readPrimitive(ifs, ngens);
        gens.resize(ngens);
        size_t n2read;
        readPrimitive(ifs, n2read);
        for (size_t i=0; i<n2read; ++i) {
            size_t igen;
            if (n2read<ngens/2) readPrimitive(ifs, igen); else igen = i;
            auto& gen = gens.at(igen);
            auto& males = gen.first;
            auto& females = gen.second;
            size_t nmales;
            size_t nfemales;
            readPrimitive(ifs, nmales);
            readPrimitive(ifs, nfemales);
            males.resize(nmales);
            females.resize(nfemales);
            for (auto* persons : { &males, &females } )
            for (auto& person : *persons) {
                readPrimitiveVector(ifs, person.chromosomeFromFather.m_intervals);
                readPrimitiveVector(ifs, person.chromosomeFromMother.m_intervals);
            }
        }
    }
}

