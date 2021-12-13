#!/usr/bin/env python
"""
syncmer class
"""

from math import isclose
import numpy as np
from jls_syncmer import Syncmer

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.
class Syncmer_Closed(Syncmer): 
    def __init__(self, k, s): # k-mer with s-minimizer at 0-offset index 0 or u
        super().__init__(k, s)
    # Returns syncmer indicator for index i.
    def indicator(self, codes):
        u = self.u
        assert len(codes) == u+1
        minimizer_value = min(codes[0],codes[u])
        for j in range( 1, u ):
            if codes[j] < minimizer_value:
                return 0
        return 1
    # Returns estimates of expectation_product_indicators as a list[0...k).
    #     E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > u=k-s.
    def expectation_product_indicator(self,i):
        u = self.u
        p = 2.0/(u+1)
        if i == 0:
            return p
        elif i < u:
            return 4.0/(i+u+1)/(u+1)
        elif i == u:
            return 1.0/(2*u+1) + 4.0/(i+u+1)/(u+1)
        return p*p
    # Returns the syncmer first passage probabilities, aka, nearest-neighbor distance distribution.
    #     P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ] 
    def first_passage_probability(self,i): # [0, 1/3, 4/21, 5/42, 5/14, 0] for k=6, u=4 
        u = self.u
        p = self.probability()
        if i == 0:
            return 0.0  
        elif i < u:
            return 4.0*u/(i+u+1)/(i+u)/(i+u-1)/p
        elif i == u:
            return (1.0/(2*u+1) + 4.0*u/(i+u+1)/(i+u)/(i+u-1))/p
        return 0.0
    # Returns the probability of an alpha-test window containing a closed syncmer.
    def to_test_probability_analytic(u,alpha):
        if alpha <= u:
            return (2*alpha)/(u+alpha)
        return 1.0

def _test_Syncmer_Closed():
    syncmer_closed = Syncmer_Closed(6,2) # Each syncmer has k-s+1 = u+1 codes.
    assert syncmer_closed.name() == 'Syncmer_Closed' 
    (syncmer_closed.get_k(),syncmer_closed.s,syncmer_closed.u) == (6,2,4)  
    (k,u) = (6,4) 
    # Tests indicator() and indicators() by hand calculation.
    assert syncmer_closed.indicator((0,1,2,3,4.0)) == 1
    assert syncmer_closed.indicator((1,2,3,4,0)) == 1
    assert syncmer_closed.indicator((1,0,2,4,3)) == 0
    assert syncmer_closed.indicators((1,2,3,4,0,5)) == [1,0]
    assert syncmer_closed.indicators((-1,1,2,3,4,0,5)) == [1,1,0]
    # Tests probability() by hand calculation.
    p0 = 0.4 # syncmer density0
    p = syncmer_closed.probability() # syncmer density
    assert isclose( p, p0 )
    # Tests expectation_product_indicators by hand calculation.
    epis0 = [2/5, 2/15, 4/35, 1/10, 1/5, 4/25] # expectation_product_indicators0
    epis = [ syncmer_closed.expectation_product_indicator(i) for i in range(k) ] # expectation_product_indicators
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests covariances without mutation, cov(Y[0],Y[i]).
    covs0 = [epi-p*p for epi in epis0] # covariances0
    covs = [syncmer_closed.covariance(i) for i in range(len(covs0))] # covariances
    assert np.allclose( covs, covs0, atol=1.0e-08 )
    # Tests covariances with mutation, cov(Y[0]*Theta_bar[0],Y[i]*Theta_bar[i]).
    theta = 0.9
    theta_bar = 1.0-theta
    mutation_alone_covariances = [ theta_bar**(i+k)-theta_bar**(2*k) for i in range(k) ] # mutation covariances cov(Theta_bar[0],Theta_bar[i]) without syncmers
    mutation_covariances0 = [ a*b for a,b in zip(epis, mutation_alone_covariances) ]
    mutation_covariances = [ syncmer_closed.mutation_covariance(i,theta) for i in range(k)]
    assert np.allclose( mutation_covariances, mutation_covariances0, rtol=1.0e-05 )
    # Tests first_passage_probabilities, f[i] by hand calculation.
    fpps0 = [0, 1/3, 4/21, 5/42, 5/14, 0] # first_passage_probabilities0
    fpps = [ syncmer_closed.first_passage_probability(i) for i in range(k) ] # first_passage_probabilities
    assert isclose(sum(fpps0), 1.0)
    assert isclose(sum(fpps), 1.0)
    assert np.allclose(fpps, fpps0)
    # Tests to_test_probability_analytic with Shaw's formula.
    length = len(fpps)
    tps0 = [ Syncmer_Closed.to_test_probability_analytic(u,alpha) for alpha in range(length) ] # test_probabilities0
    tps = Syncmer.to_test_probabilities(p,fpps) # test_probabilities
    assert np.allclose(tps,tps0)
    
def _test_Syncmer_Closed_Simulate():
    r = 100000 #(Monte Carlo realizations)
    np.random.seed(31415)

    syncmer_closed = Syncmer_Closed(6,2)
    # Tests expectation_product_indicators by simulation.
    epis_mc = syncmer_closed.expectation_product_indicators_monte_carlo(r) 
    epis0 = [ syncmer_closed.expectation_product_indicator(i) for i in range(len(epis_mc)) ]
    assert np.allclose( epis_mc, epis0, atol=1.0e-03 )
    # Tests first_passage_probabilities by simulation.
    fpps_mc = syncmer_closed.first_passage_probabilities_monte_carlo(r)
    fpps0 = [ syncmer_closed.first_passage_probability(i) for i in range(len(fpps_mc)) ]
    assert isclose(sum(fpps_mc), 1.0)
    assert isclose(sum(fpps0), 1.0)
    assert np.allclose(fpps_mc, fpps0, atol=1.0e-03 )

def main(): 
    _test_Syncmer_Closed()
    _test_Syncmer_Closed_Simulate()
    
if __name__ == "__main__":
    main()
