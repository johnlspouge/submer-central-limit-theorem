#!/usr/bin/env python
"""
syncmer class
"""

from math import isclose
from scipy.special import comb
import numpy as np
from mpmath import nsum, inf
from jls_submer import Submer

# M.C. Frith, L. Noé, G. Kucherov (2020) 
# Minimally-overlapping words for sequence similarity search. 
# Bioinformatics 36: 5344-50.

# Permits patterns of different length up to k.
# Raises an exception if the patterns overlap.
class Non_Overlapping_Pattern_Prob(Submer): 
    def __init__(self, k, length2probability): # k-mer is the submer length; length2probability = {int length (of patterns): float probability (total over patterns)}
        super().__init__(k)
        Non_Overlapping_Pattern_Prob._check(k, length2probability)
        self._length2probability = {length:p for (length,p) in length2probability.items() if p != 0.0} # Eliminates 0.0 probabilities from dictionary.
        self._prob = [] # probabilities for each pattern length
        for i in range(0,k+1):
            if i in self._length2probability:
                self._prob.append(self._length2probability[i])
            else:
                self._prob.append(0.0)
        self._cdf = [] # cumulative probabilities for each pattern length
        cdf = 0.0
        for i in range(0,k+1):
            cdf += self._prob[i]
            self._cdf.append(cdf)
        self._first_passage_probability = {}
        self._first_passage_probability_downsampled = {} # _first_passage_probability[d][i] = f^(\delta)[i]
    def get_length2probability(self):
        return self._length2probability
    def _check(k, length2probability):
        if k < max(length2probability.keys()):
            raise Exception(f'pattern longer than submer length: {k} < {max(length2probability.keys())}')
        if not isinstance(length2probability, dict):
            raise TypeError(f'length2probability is not a dictionary : {length2probability}')
        for length,probability in length2probability.items():
            if not isinstance(length, int):
                raise TypeError(f'length is not an integer in length2probability : {length}')
            if not isinstance(probability, float):
                raise TypeError(f'probability is not a float in length2probability : {probability}')
            elif probability < 0.0:
                raise Exception(f'probability < 0.0 in length2probability: {probability}')
        if 1.0 <= sum(length2probability.values()):
            raise TypeError(f'1.0 <= sum(length2probability.values()): {sum(length2probability.values())} for {length2probability.values()}')
    # Returns expectation_product_indicator E[ Y[0]*Y[i] ], e.g., E[ Y[0] ]*E[ Y[i] ] for k <= i.
    def expectation_product_indicator(self,i):
        k = self.k
        p = self._cdf[k]
        if i < 0:
            return 0.0
        elif i == 0:
            return p
        elif i < k:
            return p*self._cdf[i]
        return p*p
    # Returns the pattern first passage probabilities, aka, nearest-neighbor distance.
    #     P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ] 
    def first_passage_probability(self,i):
        k = self.k
        if i <= 0:
            return 0.0
        # Returns E[ Y[i] | Y[0]=1], the first_passage_probability.
        if i not in self._first_passage_probability:
            fpp = 0.0
            if 0 <= i <= k:
                fpp += self._prob[i]
            fpp += self.first_passage_probability(i-1)
            for n in range(1,min(i,k)+1):
                fpp -= self.first_passage_probability(i-n)*self._prob[n]
            self._first_passage_probability[i] = fpp
        return self._first_passage_probability[i]
    # Returns the pattern first passage probabilities, aka, nearest-neighbor distance.
    #     P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ] for all patterns the same length.
    def first_passage_probability_analytic(self,i):
        assert len(self._length2probability) == 1
        (m0, ph) = list(self._length2probability.items())[0]
        f = 0.0
        for j in range(0,(i-m0)//(m0-1)+1):
            f += ((-1)**(j%2))*(ph**(j+1))*comb((i-m0)-(m0-1)*j, j)
        return f
    # Returns the first_passage_probability for submers randomly down-sampled by a factor of d.
    #   d-1 = the expected number of submers before a sampled submer, consistent with Edgar's notation.
    def first_passage_probability_downsampled(self,i,d=None):
        if d is None:
            return self.first_passage_probability(i)
        if d < 1.0:
            raise Exception(f'd < 1.0 : {d}')
        if i == 0:
            return 0.0
        delta = 1.0-1.0/d
        assert 0.0 < delta < 1.0  
        fppd = self._first_passage_probability_downsampled
        # Returns the _first_passage_probability_downsampled(i,d).
        if d not in fppd:
            fppd[d] = {}
        if i not in fppd[d]:
            f_d = 0.0
            for j in range(1,i):
                f_d += self.first_passage_probability(j) * self.first_passage_probability_downsampled(i-j,d)
            f_d *= delta
            f_d += (1.0-delta)*self.first_passage_probability(i)
            self._first_passage_probability_downsampled[d][i] = f_d 
        return self._first_passage_probability_downsampled[d][i]
    # Frees the storage consumed in computing first_passage_probability_downsampled for downsampling factor d.
    def first_passage_probability_downsampled_pop(self,d):
        self._first_passage_probability_downsampled.pop(d, None)
    # Returns the first 0,1,2 moments of the submer first passage probabilities, aka, nearest-neighbor distance.
    def first_passage_moment_analytic(self,m):
        if m == 0:
            return 1.0
        elif m == 1:
            return 1.0/self.probability()
        elif m == 2:
            mu = 1.0/self.probability()
            h_prime = sum([i*self._prob[i] for i in range(len(self._prob))]) # H'(1)
            return 2.0*(1-h_prime)*mu*mu+mu
        else:
            raise Exception(f'm in [0,1,2] : m = {m}.')
    # Returns the pattern alpha-test probabilities 
    #   from generalization of analytic formula in Shaw & Yu SI p.7.
    def alpha_test_probability_analytic(self,alpha):
        assert len(self._length2probability) == 1
        (m0, ph) = list(self._length2probability.items())[0]
        n = m0-1
        pr_f_alpha = 0.0
        for i in range(1,alpha+1):
            pr_f_alpha += ((-1)**((i+1)%2))*(ph**i)*comb(alpha-(i-1)*n, i)
        return pr_f_alpha

def _test_Non_Overlapping_Pattern_Prob():
    # Checks first_passage_probability against first_passage_probability_analytic.
    non_overlapping_pattern_prob = Non_Overlapping_Pattern_Prob(3,{2:0.01})
    nopp = non_overlapping_pattern_prob
    MAX = 10
    fpps = [ nopp.first_passage_probability(i) for i in range(MAX) ] # first_passage_probabilities
    fpps0 = [0.0, 0.0, 0.01, 0.01, 0.0099, 
             0.00980000, 0.00970100, 0.009603, 0.00950599, 0.00940996]
    assert np.allclose( fpps, fpps0, atol=1.0e-08 )
    fpps0 = [ nopp.first_passage_probability_analytic(i) for i in range(MAX) ] # first_passage_probabilities
    assert np.allclose( fpps, fpps0, atol=1.0e-08 )
    # Checks conversion of first_passage_probability to alpha-test probabilities.
    tps = Submer.to_test_probabilities(nopp.probability(), fpps)
    tps0 = [0.0, 0.01, 0.02, 0.0299, 0.0397, 
            0.049401, 0.059004, 0.06851, 0.07791995, 0.0872348501]
    assert np.allclose( tps, tps0, atol=1.0e-08 )
    tps0 = [ nopp.alpha_test_probability_analytic(alpha) for alpha in range(MAX) ]
    assert np.allclose( tps, tps0, atol=1.0e-08 )
    # Checks general first_passage_probability.
    non_overlapping_pattern_prob = Non_Overlapping_Pattern_Prob(3,{1:0.1,2:0.01}) # {a,c,g,t} "a","cg"
    nopp = non_overlapping_pattern_prob
    k = 3
    assert nopp.get_k() == k
    assert np.allclose(nopp._prob, [0.0,0.1,0.01,0.0])
    assert np.allclose(nopp._cdf, [0.0,0.1,0.11,0.11])
    # Tests probability() by hand calculation.
    p0 = 0.11 # submer density0
    p = nopp.probability() # submer density
    assert isclose( p, p0 )
    # Tests expectation_product_indicators by hand calculation.
    epis0 = [0.11, 0.011, 0.0121] # expectation_product_indicators0
    epis = [ nopp.expectation_product_indicator(i) for i in range(k) ] # expectation_product_indicators
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests covariances without mutation, cov(Y[0],Y[i]).
    covs0 = [epi-p*p for epi in epis0] # covariances0
    covs = [nopp.covariance(i) for i in range(len(covs0))] # covariances
    assert np.allclose( covs, covs0, atol=1.0e-08 )
    # Tests covariances with mutation, cov(Y[0]*Theta_bar[0],Y[i]*Theta_bar[i]).
    theta = 0.9
    theta_bar = 1.0-theta
    mutation_alone_covariances = [ theta_bar**(i+k)-theta_bar**(2*k) for i in range(k) ] # mutation covariances cov(Theta_bar[0],Theta_bar[i]) without syncmers
    mutation_covariances0 = [ a*b for a,b in zip(epis, mutation_alone_covariances) ]
    mutation_covariances = [ nopp.mutation_covariance(i,theta) for i in range(k)]
    assert np.allclose( mutation_covariances, mutation_covariances0, rtol=1.0e-05 )
    # Tests first_passage_probabilities, f[i] regressively.
    MAX = 10
    fpps = [ nopp.first_passage_probability(i) for i in range(MAX) ] # first_passage_probabilities
    fpps0 = [0.0,    0.1,      0.1,       0.089,     0.0791, 
             0.0703, 0.062479, 0.0555281, 0.0493505, 0.043860169]
    assert np.allclose( fpps, fpps0, atol=1.0e-08 )
    # Tests first_passage_moment, by hand calculation and with analytic formulas.
    assert isclose(nopp.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(nopp.probability() * nopp.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    for i in range(3):
        assert isclose(nopp.first_passage_moment(i), nopp.first_passage_moment_analytic(i)) # sum(fpps0)
    assert isclose(nopp.first_passage_moment_analytic(0), 1.0)
    assert isclose(nopp.first_passage_moment_analytic(1) * nopp.probability(), 1.0)
    fpm = nsum(lambda m: m*m*nopp.first_passage_probability(m), [1, inf])
    assert isclose(nopp.first_passage_moment_analytic(2), fpm)
    # Tests first_passage_probability_downsampled.
    MAX = 20
    D = 1/0.9
    # Tests results by hand regressively.
    fppds0 = [0.0, 0.09, 0.0909, 0.081909, 0.07371908999999999, 
              0.06634719089999998, 0.05971247190899999, 0.053741224719089986, 0.04836710224719089, 0.04353039202247191]
    fppds = [ nopp.first_passage_probability_downsampled(i,D) for i in range(len(fppds0)) ]
    assert np.allclose( fppds, fppds0, atol=1.0e-08 )
    # Checks first_passage_probability_downsampled_pop, which frees memory after down-sampling computation.
    fppds = [ nopp.first_passage_probability_downsampled(i,2*D) for i in range(len(fppds0)) ]
    nopp.first_passage_probability_downsampled_pop(D)
    nopp.first_passage_probability_downsampled_pop(2*D)
    nopp.first_passage_probability_downsampled_pop(3*D) # Behaves well on unknown down-sampling factors.
    assert nopp._first_passage_probability_downsampled == {}

def main(): 
    _test_Non_Overlapping_Pattern_Prob()
    
if __name__ == "__main__":
    main()
