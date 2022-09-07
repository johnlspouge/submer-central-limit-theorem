#!/usr/bin/env python
"""
Submer class
"""

from abc import ABC, abstractmethod
from numpy import allclose
import mpmath
from math import isclose

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.

class Submer(ABC):
    def __init__(self, k):
        super().__init__()
        Submer._check(k)
        self.k = k
    def _check(k):
        if not isinstance(k, int) or k <= 0:
            raise Exception(f'not isinstance(k, int) or k <= 0: {k}')
    # Returns submer name.
    def name(self):
        return type(self).__name__
    def get_k(self):
        return self.k
    # Returns the submer probability (a.k.a. density) as E[ Y[0] ]*E[ Y[0] ].
    def probability(self):
        return self.expectation_product_indicator(0)
    @abstractmethod
    # Returns expectation_product_indicator E[ Y[0] ]*E[ Y[i] ].
    #     E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > w=k-s.
    def expectation_product_indicator(self,i):
        pass
    # Returns the submer mutation covariances.
    def mutation_covariance(self,i,theta):
        k = self.k
        theta_bar = 1.0-theta
        return self.expectation_product_indicator(i)*(theta_bar**(i+k)-theta_bar**(2*k))
    @abstractmethod
    # Returns the submer first passage probabilities, aka, nearest-neighbor distance.
    #     P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ] 
    def first_passage_probability(self,i):
        pass
    # Returns the moments of the submer first passage probabilities, aka, nearest-neighbor distance.
    def first_passage_moment(self,m,
                             k=0): # k is the shift in the island statistics
        assert isinstance(m, int) and 0 <= m
        assert isinstance(k, int) and 0 <= k
        return mpmath.nsum(lambda i: self.first_passage_probability(int(i+k))*(i**m), [1, mpmath.inf] )
    # Returns the submer mutation covariances.
    def covariance(self,i):
        p = self.probability()
        return self.expectation_product_indicator(i)-p*p

    # Jim Shaw's formula is tested by:
    #   jls_syncmer_closed._test_Syncmer_Closed_Shaw()
    #   jls_syncmer_open._test_Syncmer_Open_Shaw()

    # Returns the list of test_probabilities calculated from the list of submer first_passage_probabilities.
    #    The submer first passage probabilities are P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ].
    #    len(test_probabilities) == len(first_passage_probabilities).
    def to_test_probabilities(
        probability, # submer density
        first_passage_probabilities): # list of first-passage probabilities
        mu = probability
        fpps = first_passage_probabilities
        assert fpps[0] == 0.0
        sum0 = sum1 = 0.0
        tps = [0.0] # test_probabilities
        for alpha in range(1,len(fpps)):
            sum0 += fpps[alpha]
            sum1 += alpha*fpps[alpha]
            total = sum1+alpha*(1.0-sum0)
            tps.append(total)
        return [tp*mu for tp in tps]
    # The inversion of Jim Shaw's formula
    # Returns the submer first passage probabilities P[ Y[i]=1 and Y[j]=0 for 0<j<i | Y[0]=1 ].
    #    len(first_passage_probabilities) = len(test_probabilities) - 1.
    def to_first_passage_probabilities(
        probability, # submer density
        test_probabilities): # list of test probabilities
        mu = probability
        tps = test_probabilities
        assert tps[0] == 0.0
        fpps = [0.0] # first passage probabilities
        for alpha in range(1,len(tps)-1):
            fpps.append(Submer.to_first_passage_probability(mu,tps[alpha-1:alpha+2]))
        return fpps
    # Returns the submer first passage probability f[i] from a slice test_probabilities[i-1:i+2].
    def to_first_passage_probability(
        probability, # submer density
        test_probabilities): # slice from test_probabilities of length = 3
        mu = probability
        tps = test_probabilities
        assert len(tps) == 3
        f = -(tps[0]-2.0*tps[1]+tps[2])/mu
        return f
    # J Shaw & Y-W Yu (2021)
    # Theory of local k-mer selection with applications to long-read alignment
    # Bioinformatics DOI: 10.1093/bioinformatics/btab790
    # Theorem 4
    # Consider the simple mutation model with mutation probability per letter theta.
    # Returns the probability that 2k-1 letters contain a run of alpha=beta+k unmutated letters (0 <= beta < k).
    def prob_unmutated_run(k, theta, alpha):
        beta = alpha-k
        assert 0 <= beta < k
        if beta == k-1:
            return (1.0-theta)**(2*k-1)
        return (2.0+(k-beta-2)*theta)*theta*(1.0-theta)**(k+beta)
    # Tests prob_unmutated_run().
    def _test_prob_unmutated_run():
        k = 4
        theta = 0.1
        pur0 = [ 0.14434200000000003, 0.12400290000000003, 0.10628820000000001, 0.4782969000000001 ]
        pur = [ Submer.prob_unmutated_run(k, theta, k+beta) for beta in range(k) ]
        assert allclose(pur, pur0)
    
# Tests Shaw's formula and its inversion interconverting test_probabilities Pr(F,alpha) and first_passage_probabilities, f[i].
def _test_Syncmer_Closed_Shaw():
    u = 4 # Syncmer_Closed(6,2)
    p = 2.0/(u+1) # syncmer_closed.probability()
    # Tests conversion of first passage probabilities to test_probabilities with Shaw's formula.
    fpps0 = [0, 1/3, 4/21, 5/42, 5/14, 0] # first_passage_probabilities for Syncmer_Closed(6,2)
    tps0 = [0, 2/5, 2/3, 6/7, 1, 1] # first_passage_probabilities for Syncmer_Closed(6,2)
    tps = Submer.to_test_probabilities(p,fpps0) # test_probabilities
    assert allclose(tps, tps0) 
    # Tests conversion of test_probabilities back to first passage probabilities with Shaw's formula.
    fpps = Submer.to_first_passage_probabilities(p, tps)
    assert allclose(fpps, fpps0[0:-1]) 
    # Tests single conversion of test_probabilities back to first passage probabilities with Shaw's formula.
    fpp = Submer.to_first_passage_probability(p, tps0[1:4]) # [2/5, 2/3, 6/7]
    fpp0 = fpps0[2] # 4/21
    assert isclose(fpp,fpp0)

def main(): 
    _test_Syncmer_Closed_Shaw()
    Submer._test_prob_unmutated_run()
    
if __name__ == "__main__":
    main()
