#!/usr/bin/env python
"""
Syncmer_Open class
"""

from math import isclose
import numpy as np
from jls_syncmer import Syncmer

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.
class Syncmer_Open(Syncmer): 
    def __init__(self, k, s, t): # k-mer with s-codes whose minimum is at 0-offset index t
        self.t = t
        self._first_passage_probability = dict() # _first_passage_probability[i] = f[i]
        self._fpp = { # Stores recursion values for first_passage_probability.
            'right':dict(),
            'left':dict(),
            'zero':dict(),
        }
        super().__init__(k, s)
    # Returns syncmer indicator for index i.
    def indicator(self, codes):
        u = self.u
        t = self.t
        assert len(codes) == u+1
        minimizer = t
        minimizer_value = codes[minimizer]
        for j in range(u+1):
            if j != minimizer and codes[j] < minimizer_value:
                return 0
        return 1
    # Returns estimates of expectation_product_indicators as a list[0...k).
    #     E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > u=k-s.
    def expectation_product_indicator(self,i):
        t = self.t
        u = self.u
        p = 1.0/(u+1)
        if u-t < t:
            t = u-t
        if i == 0:
            return p
        elif i <= t:
            return 0.0
        elif i <= u-t:
            return 1.0/(i+u+1)/(u+1)
        elif i <= u:
            return 2.0/(i+u+1)/(u+1)
        return p*p
    # Returns the syncmer first passage probabilities, aka, nearest-neighbor distance.
    #     P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ] 
    def first_passage_probability(self,i):
        u = self.u
        t = self.t
        p = self.probability()
        right = self._fpp['right']
        left = self._fpp['left']
        zero = self._fpp['zero']
        if u-t < t:
            t = u-t
        # Returns the probability for only a right (k,s,t)-syncmer in the interval [0:a].
        def _right_recursion(a):
            if a < u:
                return 0.0
            elif a not in right:
                ra = _zero_recursion(a-u+t-1)
                for j in range(t):
                    ra += _right_recursion(a-j-1)
                right[a] = ra/(a+1)
            return right[a]
        # Returns the probability for only a left (k,s,t)-syncmer in the interval [0:a].
        def _left_recursion(a):
            if a < u:
                return 0.0
            elif a not in left:
                ra = _zero_recursion(a-t-1)
                for j in range(u-t):
                    ra += _left_recursion(a-j-1)
                left[a] = ra/(a+1)
            return left[a]
        # Returns the probability for no (k,s,t)-syncmer in the interval [0:a].
        def _zero_recursion(a):
            if a < u:
                return 1.0
            elif a < 0:
                return 0.0
            elif a not in zero:
                ra = 0
                for j in range(t):
                    ra += _zero_recursion(a-j-1)
                for j in range(u-t):
                    ra += _zero_recursion(a-j-1)
                zero[a] = ra/(a+1)
            return zero[a]
        # Returns E[ Y[i] | Y[0]=1], the first_passage_probability.
        if i not in self._first_passage_probability:
            fpp = (_right_recursion(i+u-t-1)+_left_recursion(i+t-1))/(i+u+1)
            self._first_passage_probability[i] = fpp/p # Divides by p for conditional probability
        return self._first_passage_probability[i]

def _test_Syncmer_Open():
    syncmer_open = Syncmer_Open(6,2,3)
    assert syncmer_open.name() == 'Syncmer_Open' 
    (syncmer_open.get_k(),syncmer_open.s,syncmer_open.t,syncmer_open.u) == (6,2,3,4) 
    k = 6
    # Tests probability() by hand calculation.
    p0 = 0.2 # syncmer density0
    p = syncmer_open.probability() # syncmer density
    assert isclose( p, p0 )
    # Tests indicator() and indicators() by hand calculation.
    assert syncmer_open.indicator((3,1,2,0,4.0)) == 1
    assert syncmer_open.indicator((1,2,3,0,4)) == 1
    assert syncmer_open.indicator((1,0,2,4,3)) == 0
    assert syncmer_open.indicators((1,2,3,0,4,5)) == [1,0]
    assert syncmer_open.indicators((-1,1,2,3,0,4,5)) == [0,1,0]
    # Tests expectation_product_indicators regressively.
    epis0 = [0.2, 0.0, 0.02857143, 0.025, 0.04444444, 0.04] # expectation_product_indicators0
    assert len(epis0) == k
    epis = [ syncmer_open.expectation_product_indicator(i) for i in range(k) ] # expectation_product_indicators
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests covariances without mutation, cov(Y[0],Y[i]).
    covs0 = [epi-p*p for epi in epis0] # covariances0
    covs = [syncmer_open.covariance(i) for i in range(len(covs0))] # covariances
    assert np.allclose( covs, covs0, atol=1.0e-08 )
    # Tests covariances with mutation, cov(Y[0]*Theta_bar[0],Y[i]*Theta_bar[i]).
    theta = 0.9
    theta_bar = 1.0-theta
    mutation_alone_covariances = [ theta_bar**(i+k)-theta_bar**(2*k) for i in range(k) ] # mutation covariances cov(Theta_bar[0],Theta_bar[i]) without syncmers
    mutation_covariances0 = [ a*b for a,b in zip(epis, mutation_alone_covariances) ]
    mutation_covariances = [ syncmer_open.mutation_covariance(i,theta) for i in range(k)]
    assert np.allclose( mutation_covariances, mutation_covariances0, rtol=1.0e-05 ) # Are the computations consistent?
    mutation_covariances0 = [2.0e-07, 0.0, 2.85686e-10, 2.49750e-11, 4.4e-12, 3.6e-13]
    assert np.allclose( mutation_covariances, mutation_covariances0, rtol=1.0e-05 ) # regressive test
    # Tests first_passage_moment, by hand calculation.
    assert isclose(syncmer_open.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(syncmer_open.probability() * syncmer_open.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    assert isclose(syncmer_open.first_passage_moment(2), 30.2555844444103) # sum((i**2)*fpps0[i])

def _test_Syncmer_Open_Simulate():
    r = 100000 #(Monte Carlo realizations)
    np.random.seed(31415)

    syncmer_open = Syncmer_Open(6,2,4)
    # Tests expectation_product_indicators by simulation.
    # Tests expectation_product_indicators by simulation.
    epis_mc = syncmer_open.expectation_product_indicators_monte_carlo(r)
    epis = [ syncmer_open.expectation_product_indicator(i) for i in range(len(epis_mc)) ] # length k
    assert np.allclose( epis_mc, epis, atol=2.0e-03 )
    u = 6-2 # k-s
    assert len(epis_mc) == u+1 # expectation_product_indicators_monte_carlo return len(list) == k.
    # Tests first_passage_probabilities by simulation.
    fpps_mc = syncmer_open.first_passage_probabilities_monte_carlo(r)
    fpps = [ syncmer_open.first_passage_probability(i) for i in range(len(fpps_mc)) ]
    assert np.allclose( fpps_mc, fpps, atol=3.0e-03 )
    assert isclose(sum(fpps_mc),1.0)
    # Tests first_passage_probabilities by simulation.
    syncmer_open = Syncmer_Open(7,2,3) # Checks that a middle position is still treated correctly.
    fpps_mc = syncmer_open.first_passage_probabilities_monte_carlo(r)
    fpps = [ syncmer_open.first_passage_probability(i) for i in range(len(fpps_mc)) ]
    assert np.allclose( fpps_mc, fpps, atol=2.0e-03 )
    assert isclose(sum(fpps_mc),1.0)

def main():
    _test_Syncmer_Open()
    _test_Syncmer_Open_Simulate()
    
if __name__ == "__main__":
    main()
