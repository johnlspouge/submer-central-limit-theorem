#!/usr/bin/env python
"""
Minimizer class.
"""

#from abc import ABC, abstractmethod
import numpy as np
from math import isclose
from jls_submer import Submer

# S. Schleimer et al. (2003)
# Winnowing: Local Algorithms for Document Fingerprinting
# SIGMOD 2003 76-85

# M.Roberts et al. (2004) 
# Reducing storage requirements for biological sequence comparison
# Bioinformatics 15:3363-3369.

class Minimizer(Submer):
    def __init__(self, w, k=1): # (w,k)-minimizer
        self.w = w
        self.k = k # Only the lower bound for genome length requires the k-mer length.
        self._store = {
            'expectation_product_indicator':dict(),
            'naught':dict(),
            'common_submer_one_sided':dict(),
        }
        super().__init__(k)
    def get_k(self):
        return super().get_k()
    # Returns a complete set of indicators.
    #   Updates the indicators by writing 1 at all syncmers after last.
    #   Assumes the last syncmer computed is codes[last].
    #   Assumes len(indicators) == len(codes).
    def indicators(self, codes):
        w = self.w
        if len(codes) < w:
            return None
        # Updates the indicators with last = index of last minimizer.
        indicators = [0]*len(codes)
        last = codes.index(min(codes[0:w])) # the index of the first minimizer in codes
        indicators[last] = 1
        last = self.update_indicators(indicators,codes,last)
        return indicators
    #   Assumes the last minimizer computed is codes[last].
    def update_indicators(self,indicators,codes,last):
        assert len(indicators) == len(codes) 
        assert last < len(indicators)
        assert indicators[last] == 1
        # Updates the indicators by writing 1 at all minimizers after last.
        #   Assumes codes[last] is the code at the last minimizer.
        def _distance_to_next_minimizer(codes):
            m = min(w+1,len(codes))
            for j in range(1,m):
                if codes[j] < codes[0]:
                    return j
            if m < w+1:
                return None
            return codes.index(min(codes[1:w+1]))
        w = self.w
        while True: # the last minimizer
            distance = _distance_to_next_minimizer(codes[last:])
            if distance is None:
                for i in range(last+1,len(codes)):
                    indicators[i] = 0
                return last # Returns the last minimizer in indicators
            for i in range(last+1,last+distance):
                indicators[i] = 0
            last += distance
            indicators[last] = 1
        
    # Returns the minimizer probability (a.k.a. density).
    def probability(self):
        return self.expectation_product_indicator(0)
    # Returns expectation_product_indicator E[ Y[0]*Y[i] ].
    #     In fact, E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i >= 2w-1.
    def expectation_product_indicator(self,i):
        # Returns the probability that Y[0]=Y[i]=1 through minimizers within (-w:i+w).
        w = self.w
        p = 2.0/(w+1)
        if i == 0:
            return p
        elif 2*w-1 <= i: # The windows [-w+1:w-1] and [i-w+1:i+w-1] share no code.
            return p*p
        # Requires calculation for i > 0.
        two_together = dict()
        naught = self._store['naught']
        # Returns the probability that Y[0]=Y[i]=1 through minimizers within (-a:b).
        def _two_together_recursion(a,b): # i and w are passed through calling scope.
            # Avoids examining a longer interval than necessary for minimizers at 0 and i.
            if a+b <= w:
                return 0.0
            # the algorithm
            if (a,b) not in two_together:
                ra = sum([_two_together_recursion(j,b) for j in range(1,a)])
                ra += _naught_recursion(i,b-i)
                ra += sum([_naught_recursion(a,j)*_naught_recursion(i-j,b-i) for j in range(1,i)])
                ra += _naught_recursion(a,i)
                ra += sum([_two_together_recursion(a,j) for j in range(i+1,b)])
                two_together[(a,b)] = ra/(a+b-1)
            return two_together[(a,b)]
        # Returns the probability that Y[0]=1 through minimizers within (-a:b).
        def _naught_recursion(a,b): # w is passed through calling scope.
            # Avoids examining a longer interval than necessary for a minimizer at 0.
            if a+b <= w:
                return 0.0
            # the algorithm
            if (a,b) not in naught:
                ra = sum([_naught_recursion(j,b) for j in range(1,a)])
                ra += 1.0
                ra += sum([_naught_recursion(a,j) for j in range(1,b)])
                naught[(a,b)] = ra/(a+b-1)
            return naught[(a,b)]
        # Returns E[ Y[i] | Y[0]=1], the first_passage_probability.
        if i not in self._store['expectation_product_indicator']:
            self._store['expectation_product_indicator'][i] = _two_together_recursion(w,i+w)
        return self._store['expectation_product_indicator'][i]
    # Returns the expected minimizer count for alpha codes.
    def expected_count_of_minimizers(w, alpha): 
        if alpha < w:
            return 0.0
        return (2*alpha-w+1)/(w+1)
    # Returns the minimizer first passage probabilities, aka, nearest-neighbor distance.
    def first_passage_probability(self, i):
        w = self.w
        if 0 < i <= w:
            return 1.0/w
        return 0.0
    # Returns the probability of an alpha-test window containing a closed syncmer.
    def to_test_probability_analytic(w,alpha):
        if alpha <= w:
            return 1.0-(w-alpha+1)*(w-alpha)/(w+1)/w
        return 1.0
    # Returns estimates of expectation_product_indicators as a list[0...length).
    #     Note E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i >= 2w-1.
    def expectation_product_indicators_monte_carlo(self, r): # r = #(realizations) 
        w = self.w
        END = 2*w
        def _update_epis(epis,indicators):
            for m in range(0,END): # epis[0] = self.probability()
                if indicators[m] == 1:
                    epis[m] += 1
        epis = [0]*END # expectation_product_indicators conditioned on initial minimizer
        WINDOW = END+w
        codes = list(np.random.uniform(size=WINDOW))
        indicators = self.indicators(codes)
        last = len(indicators)-indicators[::-1].index(1)-1 # index of last 1 in indicators
        # state has a minimizer at index 0, but it is not stationary.
        codes = codes[last:]+list(np.random.uniform(size=last))
        indicators = indicators[last:]+([0]*last)
        last = self.update_indicators(indicators,codes,0)
        # state has a minimizer at index 0, and it is stationary.
        for total_realizations in range(r): # Counts the minimizers realized within the sequence.
            _update_epis(epis, indicators)
            m_next = 1+indicators[1:].index(1)
            codes = codes[m_next:]+list(np.random.uniform(size=m_next)) # Moves window of codes.
            indicators = indicators[m_next:]+([0]*m_next) # Moves window of indicators.
            last = self.update_indicators(indicators,codes,last-m_next) # Updates the indicators.
        mu = self.probability()
        return [mu*epi/r for epi in epis]

def _test_Minimizer():
    minimizer = Minimizer(3,2) # (w,k)
    assert minimizer.name() == 'Minimizer' 
    w = minimizer.w
    assert w == 3 
    assert minimizer.get_k() == 2 
    p0 = 0.5 # minimizer density0
    p = minimizer.probability() # minimizer density
    assert isclose( p, p0 )
    # Tests indicators() by hand calculation.
    assert minimizer.indicators((0,1,2)) == [1,0,0]
    assert minimizer.indicators((1,0,2)) == [0,1,0]
    assert minimizer.indicators((1,2,0)) == [0,0,1]
    assert minimizer.indicators((0,1,2,3,4,5)) == [1,1,1,1,0,0]
    assert minimizer.indicators((1,0,2,4,3,5)) == [0,1,1,0,1,0]
    assert minimizer.indicators((1,2,0,4,5,3)) == [0,0,1,0,0,1]
    # Tests expectation_product_indicator regressively.
    epis = [minimizer.expectation_product_indicator(i) for i in range(2*w-1)]
    epis0 = [0.5, 0.16666666666666666, 0.21587301587301586, 0.26646825396825397, 0.25337301587301586]
    assert np.allclose(epis,epis0)
    EXTRA = 10 # uncorrelated for i >= 2w-1 = 5.
    for i in range(2*w-1,2*w-1+EXTRA):
        assert isclose(minimizer.expectation_product_indicator(i), p0*p0)
    # Tests first_passage_probability().
    length = 10
    fpps = [minimizer.first_passage_probability(i) for i in range(length)]
    fpps0 = [1.0/w if 0 < i <= w else 0.0 for i in range(length)]
    assert np.allclose(fpps,fpps0)
    # Tests first_passage_moment, by hand calculation.
    assert isclose(minimizer.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(minimizer.probability() * minimizer.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    assert isclose(minimizer.first_passage_moment(2), 4.66666666666667) # sum((i**2)*fpps0[i])
    # Tests to_test_probability_analytic() with Shaw's formula.
    tps = [ Minimizer.to_test_probability_analytic(w,alpha) for alpha in range(length+1) ]
    tps0 = [0.0, 0.5, 0.8333333333333334, 1.0]
    assert np.allclose(tps[:4],tps0)
    fpps = Submer.to_first_passage_probabilities(p,tps)
    assert np.allclose(fpps,fpps0)

def _test_Minimizer_Simulate():
    r = 10000 #(Monte Carlo realizations)
    np.random.seed(31415)

    minimizer = Minimizer(3,2)
    w = minimizer.w
    LENGTH = 2*w
    epis = [minimizer.expectation_product_indicator(i) for i in range(LENGTH)]
    epis_mc = minimizer.expectation_product_indicators_monte_carlo(r)
    assert np.allclose(epis,epis_mc,atol=3.0e-03)

def main(): 
    _test_Minimizer()
    _test_Minimizer_Simulate()
    
if __name__ == "__main__":
    main()
    