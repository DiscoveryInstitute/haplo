#include "ParentPairUrn.h"
#include <random>
#include <utility>
#include "PolyaUrn.h"

/**
 *  This method chooses mothers and fathers for children according to the scheme given in 
 *  http://bio-complexity.org/ojs/index.php/main/article/view/BIO-C.2016.4 eqn. 68.
 *
 *  The variables have different names for readability in C++, but I have attempted to identify the
 *  correspondence where possible. It uses a slightly more complex approach than described for 
 *  choosing fathers (in order to be more space-efficient) and so some of the variables do not map 
 *  one-to-one.
 *
 *  const long nMaleAdults;      <-  M_{t+1} in paper
 *  const long nFemaleAdults;    <-  F_{t+1} in paper
 *  const double alpha;          <-  alpha in paper
 *  const double beta;           <-  beta in paper
 *  
 *  PolyaUrn femalesPolyaUrn;    
 *    ^ Used to choose mothers for children. Records how many children each female has already had,
 *        because these may influence probability of further children.
 *        The urn returns an id of a mother, from 0 to f_{t+1}, 
 *        where f_{t+1} is the previous number of mothers.
 *      NB: Variables f_{t+1} and C_{t+1},f from the paper are hidden in femalesPolyaUrn.
 *
 *  std::vector<Mother> mothers;
 *    ^ Represents information about mothers once they have been chosen.
 *
 *  PolyaUrn Mother.partnersPolyaUrn;    
 *    ^ Used to choose fathers for children. For each female, records how many children she has 
 *        already had with each male, as may influence probability of further children with him.
 *        The urn returns an id from 0 to m(f)_{t+1} where m(f) is the previous number of 'husbands' 
 *        this female has had, indicating whether he is her first, second, third, ...
 
 *  std::vector<long> Mother.partnerIds;
 *    ^ For each female, identify who the first, second, third, ..., 'husbands' are.
 *      NB: Variables C_{t+1},mf are implictly hidden in malesPolyaUrnForFemale and malesForFemale.
 *
 *   long nFathers=0;  <- m_{t+1} in paper
 *    ^ Records how many males have already become fathers. This is also the id for the next father.
 *
 */
using namespace std;

uniform_real_distribution<double> ParentPairUrn::random(0.0, 1.0);

ParentPairUrn::ParentPairUrn(long nM, long nF, double a, double b) :
        m_nMaleAdults(nM), m_nFemaleAdults(nF), m_alpha(a), m_beta(b),
        m_femalesPolyaUrn(m_nFemaleAdults, m_alpha, true),
        m_mothers(),
        m_nFathers(0) {
    m_mothers.reserve(nF);
}

pair<long,long> ParentPairUrn::pick(mt19937_64& rng) {

    // First choose a female mother. 
    //   Probability depends on how many children she has had (PolyaUrn remembers that).
    size_t f = m_femalesPolyaUrn.sample(rng);
    
    // If this is a new mother having the newest id, make some space for her data.
    if (f == m_mothers.size()) {
        m_mothers.push_back(Mother { PolyaUrn(m_nMaleAdults, m_beta, false), vector<long>{} });
    }
    
    // Second, choose a male father. 
    //   Probability depends on who mother is and how many children they have had together.
    //   First we choose a proxy: 0 means her first 'husband', 1 means her second 'husband' ...
    size_t mp = m_mothers[f].partnersUrn.sample(rng);
    size_t m;
    // Then convert the proxy into a reference to a particular male.
    if (mp == m_mothers[f].partnerIds.size()) { // If he didn't have a child with her before.
        long mc = (long)(random(rng) * m_nMaleAdults); // Choose who he via a random id.
        m = (mc < m_nFathers) // If this id has been used before, the male already has 
            ? mc              // children. If it has not been used before, this is a new Dad,
            : m_nFathers++;   // so give him the next unused id and record we have extra Dad.
        m_mothers[f].partnerIds.push_back(m); // Record he was the latest father with this female.
    } else { // If he did have a child with her before.
        m = m_mothers[f].partnerIds[mp]; // Just look up his id. 
        // NB: The PolyaUrn now remembers they had more than one child, but we don't need to.
    }
    
    return pair<long, long>(m,f);

}
       
            
            
            
            
    