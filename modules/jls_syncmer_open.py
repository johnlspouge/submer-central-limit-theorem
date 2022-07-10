#!/usr/bin/env python
"""
Syncmer_Open class
"""
from math import isclose
import numpy as np
from jls_syncmer_parametrized import Syncmer_Parametrized

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.
class Syncmer_Open(Syncmer_Parametrized): 
    def __init__(self, k, s, t, eps=0.0): # k-mer with s-codes whose minimum is at 0-offset index t
        ts = [t]
        self.t = t
        super().__init__(k, s, ts, eps)

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
    # Tests lagging moments for island calculations.
    syncmer_open = Syncmer_Open(6,2,4)
    k=6
    p = syncmer_open.probability() # syncmer density
    assert isclose(syncmer_open.first_passage_moment(0,k)*p, 0.0575757575757576) 
    assert isclose(syncmer_open.first_passage_moment(1,k)*p, 0.174365079365079)
    assert isclose(syncmer_open.first_passage_moment(2,k)*p, 0.800723325467849) 
    # f = [0.0, 0.16666666666666666, 0.11904761904761904, 0.08928571428571427, 0.06944444444444443, 0.15555555555555553, 0.11212121212121211, 0.0835137085137085, 0.06379731379731379, 0.04972309436595151, 0.03272297808012093, 0.021625030062530057, 0.014168337789661319, 0.009049782824292626, 0.0054801691393170085, 0.003296003321379261, 0.0019523349068417092, 0.0011308345612351303, 0.0006390323026802223, 0.00035743922998532834, 0.00019700395873133167, 0.00010676000395039424, 5.696909398239726e-05, 3.0077290745075236e-05]
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
