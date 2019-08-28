#include "PopulationStructureHistory.h"
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

using namespace std;

PopulationStructureHistory PopulationStructureHistory::create(
        string type, long initial_gen, long init_popn, long final_popn, double growth_rate) {
    if (type == "Constant") return Constant(initial_gen, final_popn);
    if (type == "Linear") return Linear(initial_gen, init_popn, final_popn);
    if (type == "Exponential") return Exponential(initial_gen, init_popn, final_popn);
    if (type == "Logistic") return Logistic(initial_gen, init_popn, final_popn, growth_rate);
    if (type == "Sinusoidal") return Sinusoidal(initial_gen, init_popn, final_popn, growth_rate);
    throw invalid_argument("Unknown type argument to PSH");
}

inline pair<long, long> evenSplit(long popn) { return pair<long, long> { popn/2, popn/2 + popn%2}; }

inline long nearestLong(double x) { return (long)(x>0 ? x+0.5 : x-0.5); }

PopulationStructureHistory PopulationStructureHistory::Constant(
        long initial_gen, long males, long females) {
	pair<long, long> nMaleFemale(males, females);
    return PopulationStructureHistory(vector<pair<long,long>>(initial_gen+1, nMaleFemale));
}

PopulationStructureHistory PopulationStructureHistory::Constant(long initial_gen, long popn) {
    return PopulationStructureHistory(vector<pair<long,long>>(initial_gen+1, evenSplit(popn)));
}

PopulationStructureHistory PopulationStructureHistory::Linear(
        long initial_gen, long initial_popn, long final_popn) {
	vector<pair<long,long>> vec;
    vec.reserve(initial_gen+1);
    double step = (initial_popn - final_popn) / (double) initial_gen;
    for(int i=0; i<initial_gen; ++i) {
        long popn = nearestLong(final_popn + i*step);
        vec.emplace_back(evenSplit(popn));
    }
    vec.emplace_back(evenSplit(initial_popn));
	return PopulationStructureHistory(vec);
}

PopulationStructureHistory PopulationStructureHistory::Exponential(
        long initial_gen, long initial_popn, long final_popn) {
	vector<pair<long,long>> vec;
    vec.reserve(initial_gen+1);
    double popn = final_popn;
    double factor = exp( (log(initial_popn) - log(final_popn)) / initial_gen );
    for(int i=0; i<initial_gen; ++i) {
        vec.emplace_back(evenSplit(nearestLong(popn)));
		popn *= factor;
    }
    vec.emplace_back(evenSplit(initial_popn));
	return PopulationStructureHistory(vec);
}

PopulationStructureHistory PopulationStructureHistory::Logistic(
        long initial_gen, long initial_popn, long final_popn, double growth_rate) {
	vector<pair<long,long>> vec;
    vec.reserve(initial_gen+1);
    double popn = final_popn;
    for(int i=0; i<initial_gen; ++i) {
        vec.emplace_back(evenSplit(nearestLong(popn)));
        popn += growth_rate * popn * (initial_popn - popn) / initial_popn;
        if (popn > final_popn) popn = final_popn;
    }
    vec.emplace_back(evenSplit(initial_popn));
	return PopulationStructureHistory(vec);
}

PopulationStructureHistory PopulationStructureHistory::Sinusoidal(
        long initial_gen, long initial_popn, long final_popn, double growth_rate) {
	vector<pair<long,long>> vec;
    vec.reserve(initial_gen+1);
    constexpr double pi = 3.14159265;
    double base = final_popn;
    double amplitude = 0.5 * (initial_popn - final_popn);
    double k = 2 * pi * growth_rate / initial_gen;
    for(int i=0; i<initial_gen; ++i) {
		long popn = nearestLong(base + amplitude * (1 - cos(k*i)));
        vec.emplace_back(evenSplit(popn));
    }
	return PopulationStructureHistory(vec);
}

long PopulationStructureHistory::getTotalNumberOfMales() const {
    long sum = 0L;
    for (auto& popn : m_populations) sum += popn.first;
    return sum;
}

long PopulationStructureHistory::getTotalNumberOfFemales() const {
    long sum = 0L;
    for (auto& popn : m_populations) sum += popn.second;
    return sum;
}

bool PopulationStructureHistory::operator==(const PopulationStructureHistory& o) const {
    return m_populations == o.m_populations;
}

