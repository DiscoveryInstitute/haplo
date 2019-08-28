#include "PolyaUrn.h"
#include <random>
#include <vector>
#include <cstdio>
#include <cstdlib>

/** 
 *  Implements a PolyaUrn with population 'total_options' and randomness parameter 'alpha'. 
 *
 *  Previously chosen picks have additional weight proportional to the number of times they were 
 *  previously chosen. Potentially there could be a very large number of these picks, so the 
 *  cumulative weights are stored in a binary tree for easy lookup. 
 *
 *  At each level 'ilevel' of the binary tree, the relevant node is indexed by index>>ilevel
 *  (all but the bottom ilevel bits). Each node sums 2 nodes from the previous level, 
 *  or 2^ilevel nodes from the bottom layer.
 */

using namespace std;

uniform_real_distribution<double> PolyaUrn::random(0.0, 1.0); 
 
long PolyaUrn::sample(mt19937_64& rng) {
    const double alpha = m_alpha;
    const long n = m_total_options;
    const long m = m_current_options;
    const long c = m_current_picks;
    double r = random(rng);
    if (alpha <= 0.0) return 0L;
    const double nalpha = n * alpha; // <- all 'n' options have weight 'alpha'.
    r *= (nalpha + c);              // <- previously chosen options 'c' have additional unit weight.
    if (r < nalpha) {  // First part of range represents completely-random alpha-weighted choice.
        long i = (long) (r / alpha);  // Using that subrange, choose a member of the population.
        if (i >= m) return pick_new_index();           // Id has not been chosen before.
        else        return pick_by_existing_index(i);  // Id has been chosen before.
    } else {           // Second part of range represents additional weight due to previous choices.
        r -= nalpha;                  // Using that subrange, 
        long ic = (long) r;           // 'ic' represents a particular previous choice.
        return pick_by_ic(ic);        // Lookup which id it corresponds to.
    }
}


// Picking a new index is simple. The complicated part is correctly and efficiently initializing the  
// new nodes on the binary lookup tree. 
long PolyaUrn::pick_new_index() { 
    const long index = m_current_options;
    auto& nodes = m_nodes;
    int nlevels = nodes.size();
    
    ++m_current_options;
    ++m_current_picks;
    
    // Special case: initialize the first level with a single node with weight 1.
    if (nlevels == 0) {
        nodes.push_back(vector<long>{1L});
        if (m_reserve) nodes[0].reserve(m_total_options);
        return index;
    }
    
    // Special case: initialize a new level because the index has reached a new power of 2.
    // The new node carries the sum from the previous level.
    if (index == 1<<(nlevels-1)) {
        long sum = nodes[nlevels-1][0];
        nodes.push_back(vector<long>{sum});
        if (m_reserve) nodes[nlevels].reserve((m_total_options>>nlevels)+1);
        ++nlevels;
    }

    // Starting at the bottom, add any new nodes that are needed on existing levels.
    size_t i = index;
    int ilevel = 0;
    do {
        nodes[ilevel].push_back(1L);
        i>>=1;
        ++ilevel;
    } while (i == nodes[ilevel].size());
    
    // Then increment the weights on existing nodes.
    for(; ilevel<nlevels; ++ilevel) {
        ++nodes[ilevel][i]; 
        i>>=1;
    }
    return index;
}

// Picking an existing index is trivial. The slightly complicated part is updating the nodes. 
long PolyaUrn::pick_by_existing_index(const long index) {
    auto& nodes = m_nodes;
    int nlevels = nodes.size();
    
    ++m_current_picks;
    
    // This is not a new index so all the nodes already exist, need only be incremented.
    long i = index;
    for(int ilevel=0; ilevel<nlevels; ++ilevel) {
        ++nodes[ilevel][i];
        i>>=1;
    }
    return index;
}

// Looking up a previous choice is more complex. The nodes must also be updated simultaneously.
long PolyaUrn::pick_by_ic(long ic) {
    auto& nodes = m_nodes;
    int nlevels = nodes.size();
    
    ++m_current_picks;
    
    // Look up the choice in the binary tree. Start at the top and work down.
    long i = 0L;
    for (int ilevel=nlevels-1; ilevel>=0; --ilevel) {
        i <<= 1;
        if (ic >= nodes[ilevel][i]) {
            ic -= nodes[ilevel][i];
            i += 1;
        } 
        ++nodes[ilevel][i];
    }
    return i;
}
    