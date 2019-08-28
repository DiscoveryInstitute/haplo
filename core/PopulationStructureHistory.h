#pragma once
#include <string>
#include <utility>
#include <vector>

/**
    PopulationStructureHistory:
    Stores the population structure over history.
    Currently stored as a simple list of pairs: the number of males and number of females at each
    generation. Generation 0 is the extant generation from which we take the sample, last in time.
    Other generations are listed in order backwards in time.
*/
class PopulationStructureHistory {

    friend class FileCheckpointer;
    friend class FileReader;

  private:

    std::vector<std::pair<long,long>> m_populations;

	PopulationStructureHistory(std::vector<std::pair<long,long>> popns) : m_populations(popns) {}

public:

    PopulationStructureHistory() // for creating unassigned variables
        : m_populations() {}

    static PopulationStructureHistory Constant(
            long initial_gen, long final_males, long final_females);
    static PopulationStructureHistory Constant(long initial_gen, long final_popn);
    static PopulationStructureHistory Linear(
            long initial_gen, long initial_popn, long final_popn);
    static PopulationStructureHistory Exponential(
            long initial_gen, long initial_popn, long final_popn);
	static PopulationStructureHistory Logistic(
            long initial_gen, long initial_popn, long final_popn, double growth_rate);
	static PopulationStructureHistory Sinusoidal(
            long initial_gen, long initial_popn, long final_popn, double growth_rate);

	static PopulationStructureHistory create(
        std::string type, long initial_gen, long initial_popn, long final_popn, double growth_rate);

    long getNumberOfGenerations() const { return m_populations.size(); }
    long getNumberOfMales(long gen) const { return m_populations[gen].first; }
    long getNumberOfFemales(long gen) const { return m_populations[gen].second; }

    long getTotalNumberOfMales() const;
    long getTotalNumberOfFemales() const;

    bool operator==(const PopulationStructureHistory&) const;

};

