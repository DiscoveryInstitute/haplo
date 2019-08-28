#pragma once

#include <random>
#include <vector>

class PolyaUrn {
    // Parameters (see notes on constructor)
    const long m_total_options;
    const double m_alpha;
    const bool m_reserve;
    
    // State
    long m_current_options;
    long m_current_picks;
    std::vector<std::vector<long>> m_nodes;
    
    static std::uniform_real_distribution<double> random;
    long pick_by_existing_index(const long index);
    long pick_new_index();
    long pick_by_ic(long ic);
    
 public:
 
    /** 
     *  Create a new Polya Urn where the user can choose from a population of size 'total_options', 
     *  with randomness parameter 'alpha' (larger -> more random, smaller -> more concentrated).
     *  The optional parameter 'reserve' can be used to preallocate memory when a large memory 
     *  need is anticipated.
     */
    PolyaUrn(long total_options, double alpha, bool reserve = false) 
        : m_total_options(total_options), m_alpha(alpha), m_reserve(reserve),
        m_current_options(0L), m_current_picks(0L) {}
    
    /** 
     * Choose from the population using the supplied random number generator. An id is returned. 
     * The first pick is always labelled id=0. Subsequent picks are given monotonically-increasing 
     * ids if they have not been picked before. 
     */
    long sample(std::mt19937_64& rng); 
    
};
