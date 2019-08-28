#include <memory>
#include <string>
#include <utility>
#include <vector>
#include "AncestralChromosomes.h"
#include "HaploBlockVariants.h"

using namespace std;


AncestralChromosomes& AncestralChromosomes::initialiseFoundingGeneration() {
    // Initialize founding generation. Each chromosome gets its own ref_id
    const long nGenerations = m_graph->getNumberOfGenerations();
    const long foundingGeneration = nGenerations-1;

    const long gen = foundingGeneration;
    const auto& graphMaleFounders = m_graph->getMales(gen);
    const auto& graphFemaleFounders = m_graph->getFemales(gen);
    const long nMaleFounders = graphMaleFounders.size();
    const long nFemaleFounders = graphFemaleFounders.size();
    long ref = 0;
    vector<PersonNode> maleFounders;   maleFounders.reserve(nMaleFounders);
    vector<PersonNode> femaleFounders; femaleFounders.reserve(nFemaleFounders);
    for (auto nodesPair : { make_pair(&graphMaleFounders,   &maleFounders),
                            make_pair(&graphFemaleFounders, &femaleFounders) })
    for (auto& graphFounder : *nodesPair.first) {
        auto& founders = *nodesPair.second;

        auto chromosomeA = m_variants->getReferenceChromosome(ref++);
        auto chromosomeB = m_variants->getReferenceChromosome(ref++);
        chromosomeA = chromosomeA.filtered(graphFounder.chromosomeFromFatherExtancyPattern);
        chromosomeB = chromosomeB.filtered(graphFounder.chromosomeFromMotherExtancyPattern);
        founders.emplace_back(move(chromosomeA), move(chromosomeB));
    }
    m_nodes.at(gen) = { move(maleFounders), move(femaleFounders) };
    return *this;
}

AncestralChromosomes& AncestralChromosomes::forwardsDropAllGenerations(
        const double mu, mt19937_64& rng) {

    const long nGenerations = m_graph->getNumberOfGenerations();
    const long extantGeneration = 0L;
    const long foundingGeneration = nGenerations-1;

    // Drop the genes to following generations.
    for (long gen=foundingGeneration-1; gen >= extantGeneration; --gen) {
        forwardsDropOneGeneration(gen, mu, rng);
    }

    return *this;
}

AncestralChromosomes& AncestralChromosomes::forwardsDropOneGeneration(
        const long gen,
        const double mu, mt19937_64& rng) {

    const long nGenerations = m_graph->getNumberOfGenerations();
    const long lastGeneration = nGenerations-1;
    if (gen < 0L) throw invalid_argument("Gen must be >=0");
    if (gen >= lastGeneration) throw invalid_argument("Gen id too big");

    auto& maleParents = m_nodes.at(gen+1).first;
    auto& femaleParents = m_nodes.at(gen+1).second;
    const auto& graphMaleChildren = m_graph->getMales(gen);
    const auto& graphFemaleChildren = m_graph->getFemales(gen);
    long nMaleChildren = graphMaleChildren.size();
    long nFemaleChildren = graphFemaleChildren.size();
    vector<PersonNode> maleChildren;   maleChildren.reserve(nMaleChildren);
    vector<PersonNode> femaleChildren; femaleChildren.reserve(nFemaleChildren);

    for (auto nodesPair : { make_pair(&graphMaleChildren,   &maleChildren),
                            make_pair(&graphFemaleChildren, &femaleChildren) })
    for (auto& graphChild : *nodesPair.first) {
        auto& children = *nodesPair.second;

        Chromosome chromosomeFromFather;
        Chromosome chromosomeFromMother;

        if (graphChild.extantFather()) {
            auto& father = maleParents.at(graphChild.fatherId);
            chromosomeFromFather =
                    Chromosome::combine(father.chromosomeFromFather,
                                        father.chromosomeFromMother,
                                        graphChild.chromosomeFromFatherRecombinationPattern);
            chromosomeFromFather =
                    chromosomeFromFather.filtered(graphChild.chromosomeFromFatherExtancyPattern);
            chromosomeFromFather = m_variants->mutate(chromosomeFromFather, mu, gen, rng);
        }

        if (graphChild.extantMother()) {
            auto& mother = femaleParents.at(graphChild.motherId);
            chromosomeFromMother =
                    Chromosome::combine(mother.chromosomeFromFather,
                                        mother.chromosomeFromMother,
                                        graphChild.chromosomeFromMotherRecombinationPattern);
            chromosomeFromMother =
                    chromosomeFromMother.filtered(graphChild.chromosomeFromMotherExtancyPattern);
            chromosomeFromMother = m_variants->mutate(chromosomeFromMother, mu, gen, rng);
        }

        children.emplace_back(move(chromosomeFromFather), move(chromosomeFromMother));
    }
    m_nodes.at(gen) = { move(maleChildren), move(femaleChildren) };

    return *this;
}

size_t AncestralChromosomes::getSizeInMemory(long gen) const {
    size_t sum = sizeof(m_nodes.at(gen));
    for (auto* persons : { &m_nodes.at(gen).first, &m_nodes.at(gen).second })
    for (auto& person : *persons)
    for (auto* chromosome : { &person.chromosomeFromFather, &person.chromosomeFromMother }) {
        sum += chromosome->getSizeInMemory();
    }
    return sum;
}

string AncestralChromosomes::toString(long gen) const {
    std::string s;
    s = "Generation: " + to_string(gen);
    s += "\nMales: " + to_string(getMales(gen).size()) + "\n";
    for (const auto& m : getMales(gen)) {
        s += m.chromosomeFromFather.toString();
        s += m.chromosomeFromMother.toString();
    }
    s += "\nFemales: " + to_string(getFemales(gen).size()) + "\n";
    for (const auto& f : getFemales(gen)) {
        s += f.chromosomeFromFather.toString();
        s += f.chromosomeFromMother.toString();
    }
    s+="\n";
    return s;
}

bool AncestralChromosomes::operator==(const AncestralChromosomes& o) const {
    return ( m_graph == o.m_graph
             || ( m_graph != nullptr
                && o.m_graph != nullptr
                && *m_graph == *o.m_graph ) )
        && ( m_variants == o.m_variants
             || ( m_variants != nullptr
                && o.m_variants != nullptr
                && *m_variants == *o.m_variants ) )
        && m_nodes == o.m_nodes;
}

