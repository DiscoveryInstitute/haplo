#pragma once

#include <random>
#include <utility>
#include <vector>
#include "PolyaUrn.h"

// Class to choose mothers and fathers for children according to the scheme given in 
// http://bio-complexity.org/ojs/index.php/main/article/view/BIO-C.2016.4 eqn. 68.
// See more notes in cpp file.

class ParentPairUrn {

  private:
    const long m_nMaleAdults;
    const long m_nFemaleAdults;
    const double m_alpha;
    const double m_beta;
    
    struct Mother {
        PolyaUrn partnersUrn;
        std::vector<long> partnerIds;
    };
    
    PolyaUrn m_femalesPolyaUrn; 
    std::vector<Mother> m_mothers;
    long m_nFathers=0;

    static std::uniform_real_distribution<double> random;
    
  public:
    /** 
     *  Create a new ParentPairUrn containing 'nM' males (M_{t+1} in paper) 
     *  and 'nF' females (F_{t+1} in paper), with mother randomness parameter 'a' (alpha) 
     *  and father-per-mother randomness parameter 'b' (beta). 
     *  Large alpha makes child more likely to have a random mother. 
     *  Small alpha makes child more likely to have a mother who already has children.
     *  Large beta makes child more likely to have a random father.
     *  Small beta makes child more likely to have father who already has children with the mother.     
     */
    ParentPairUrn(long nM, long nF, double a, double b);
    
   /** 
     * Choose a parent-pair using the supplied random number generator. Two ids is returned, first 
     * for the father and second for the mother. The first pick is always labelled with ids=0,0. 
     * Subsequent fathers are given monotonically-increasing ids if they have not been picked 
     * before, and subsequent mothers are given monotonically-increasing ids from a separate 
     * sequence.
     */
    std::pair<long,long> pick(std::mt19937_64& rng);
    
    long getNumberOfFathers() { return m_nFathers; }
    long getNumberOfMothers() { return m_mothers.size(); }
    
};
       
            
            
            
            
    