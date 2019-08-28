#pragma once
#include <string>
#include <utility>
#include <vector>

/**
    MutationRateHistory:
    Stores the mutation rate over history.
    Currently stored as a simple list of values: the mutation rate at each generation.
    Generation 0 is the extant generation from which we take the sample, last in time.
    Other generations are listed in order backwards in time.
*/
class MutationRateHistory {

    friend class FileCheckpointer;
    friend class FileReader;

  private:

    std::vector<double> m_rates;

	MutationRateHistory(std::vector<double> rates) : m_rates(rates) {}

public:

    MutationRateHistory() {}// for creating unassigned variables

    static MutationRateHistory Constant(long ngens, double rate) {
        return MutationRateHistory(std::vector<double>(ngens,rate));
    }

    long getNumberOfGenerations() const { return m_rates.size(); }
    double getMutationRate(long gen) const { return m_rates[gen]; }

    bool operator==(const MutationRateHistory& o) { return m_rates == o.m_rates; }

};

