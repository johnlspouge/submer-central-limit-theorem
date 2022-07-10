#!/usr/bin/env python
"""
Syncmer_Open class
"""
import numpy as np
from math import isclose
from collections import Counter
from jls_syncmer import Syncmer

# A. Dutta et al (2022) 
# Parameterized syncmer schemes improve long-read mapping. 
# bioRxiv 2022.01.10.475696.

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.

# k-mer with s-codes whose minimum is at 0-offset index t from ascending tuple or list ts with rejection probability eps
class Syncmer_Parametrized(Syncmer): 
    def __init__(self,  
        k, # k-mer length
        s, # s-mer length
        ts, # s-minimizer offsets (starting at 0) that turn the k-mer into a syncmer
        eps=0.0): # rejection probability for downsampling
        self.ts = sorted(ts)
        counter = dict(Counter(self.ts))
        for e,v in counter.items():
            assert 0 <= e <= k-s
            assert v == 1 # Every position is unique.
        self.eps = eps
        # Storage of _complementary_alpha_run_probability depends on self.eps.
        self._complementary_alpha_run_probability = dict() # _complementary_alpha_run_probability[i]
        # Storage of _expectation_product_indicator does not depend on self.eps.
        self._expectation_product_indicator = dict() # expectation_product_indicator[i] = E(Y[i]|Y[0])
        super().__init__(k, s)
    def set_eps(self, eps):
        assert 0.0 <= eps < 1.0
        self.eps = eps
        self._complementary_alpha_run_probability = dict()
    # Returns indicator, without consideration of eps-downsampling.
    def indicator(self, codes):
        u = self.u
        ts = self.ts
        assert len(codes) == u+1
        min_code = min(codes)
        for t in ts:
            if codes[t] == min_code:
                return 1
        return 0
    # Returns and stores expectation_product_indicators E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > u=k-s.
    def expectation_product_indicator(self,i):
        ts = self.ts
        u = self.u
        if i < 0:
            i = -i
        assert 0 <= i
        if u < i: # The randomized s-codes are independent in [0:u] and [i:i+u].
            return self.expectation_product_indicator(0)**2
        if i not in self._expectation_product_indicator: 
            if i == 0:
                self._expectation_product_indicator[0] = len(ts)/(u+1)
            elif 0 < i <= u:
                set_ts = set(ts)
                set_ts_plus_i = set([t+i for t in ts])
                card_ts_x_ts_plus_i = len(set_ts.intersection(set_ts_plus_i))
                card_0_i_x_ts = len([ t for t in ts if t < i ]) 
                card_1u_iu_x_ts = len([ t for t in ts if u-(i-1) <= t <= u ]) 
                epi = 0.0 
                epi += card_ts_x_ts_plus_i*(u+1) 
                epi += (card_0_i_x_ts+card_1u_iu_x_ts)*len(ts)
                epi /= (i+u+1)*(u+1)
                self._expectation_product_indicator[i] = epi
        if i == 0:
            return self._expectation_product_indicator[0]*(1-self.eps)
        return self._expectation_product_indicator[i]*(1-self.eps)**2
    # Returns the syncmer complementary alpha-run probabilities.
    def complementary_alpha_run_probability(self,alpha):
        if alpha <= 0:
            return 1.0
        if alpha not in self._complementary_alpha_run_probability:
            ts = self.ts
            u = self.u
            q = 0.0
            for beta in range(u+alpha):
                c = len([t for t in ts if 0 <= beta-t < alpha])
                q += self.eps**c*self.complementary_alpha_run_probability(alpha-beta-1)*self.complementary_alpha_run_probability(beta-u)
            q /= u+alpha
            self._complementary_alpha_run_probability[alpha] = q
        return self._complementary_alpha_run_probability[alpha]
    # Returns the syncmer complementary alpha-run probabilities.
    def alpha_run_probability(self,alpha):
        return 1.0-self.complementary_alpha_run_probability(alpha)
    # Returns the syncmer complementary alpha-run probabilities.
    def first_passage_probability(self,alpha):
        if alpha == 0:
            return 0.0
        fpp = 0.0
        fpp += self.complementary_alpha_run_probability(alpha+1)
        fpp += -2.0*self.complementary_alpha_run_probability(alpha)
        fpp += self.complementary_alpha_run_probability(alpha-1)
        fpp /= self.probability()
        return fpp

def main():
    # Tests parametrized syncmers (closed syncmers) without downsampling.
    # Tests expectation_product_indicator.  .
    ts = [0,4]
    scp = Syncmer_Parametrized(6, 2, ts) # closed syncmer
    epis0 = [2/5, 2/15, 4/35, 1/10, 1/5, 4/25] # expectation_product_indicators0
    epis = [scp.expectation_product_indicator(i) for i,e in enumerate(epis0)]
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    eps = 0.9
    scp.set_eps(eps)
    epis0 = [epi0*(1.0-eps)**2 for epi0 in epis0]
    epis0[0] /= (1.0-eps) # EY[0] is the density
    epis = [scp.expectation_product_indicator(i) for i,e in enumerate(epis0)]
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests first_passage_probability.
    scp.set_eps(0.0)
    fpps0 = [0, 1/3, 4/21, 5/42, 5/14, 0] # first_passage_probabilities0
    fpps = [scp.first_passage_probability(i) for i,e in enumerate(epis0)] # first_passage_probabilities
    assert np.allclose(fpps, fpps0)
    assert isclose(sum(fpps), 1.0)
    assert isclose(scp.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(scp.probability() * scp.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    assert isclose(scp.first_passage_moment(2), 7.88095238095238) # sum((i**2)*fpps0[i])
    # Tests setting eps.
    eps = 0.9
    scp.set_eps(eps)
    fpps = [scp.first_passage_probability(i) for i in range(1000)] # first_passage_probabilities
    assert isclose(scp.first_passage_probability(1), 1/3*(1.0-eps))
    assert isclose(scp.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(scp.probability() * scp.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    # Tests parametrized syncmers (open syncmers) without downsampling.
    # Tests expectation_product_indicator.  .
    ts = [3]
    sop = Syncmer_Parametrized(6,2,ts)
    epis0 = [0.2, 0.0, 0.02857143, 0.025, 0.04444444, 0.04] # expectation_product_indicators0
    epis = [ sop.expectation_product_indicator(i) for i,e in enumerate(epis0) ] # expectation_product_indicators
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests setting eps.
    eps = 0.9
    sop.set_eps(eps)
    epis0 = [epi0*(1.0-eps)**2 for epi0 in epis0]
    epis0[0] /= (1.0-eps) # EY[0] is the density
    epis = [sop.expectation_product_indicator(i) for i,e in enumerate(epis0)]
    assert np.allclose( epis, epis0, atol=1.0e-08 )
    # Tests first_passage_probability.
    sop = Syncmer_Parametrized(6,2,ts)
    assert isclose(sop.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(sop.probability() * sop.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    assert isclose(sop.first_passage_moment(2), 30.2555844444103) # sum((i**2)*fpps0[i])
    # Tests first_passage_probability with eps-downsampling.
    eps=0.9
    sop9 = Syncmer_Parametrized(6,2,ts,eps)
    assert isclose(sop9.first_passage_probability(1), sop.first_passage_probability(1)*(1.0-eps))
    assert isclose(sop9.first_passage_moment(0), 1.0) # sum(fpps0)
    assert isclose(sop9.probability() * sop9.first_passage_moment(1), 1.0) # sum(i*fpps0[i])
    # Tests first_passage_probability.
    ts = [4]
    sop = Syncmer_Parametrized(6,2,ts)
    assert isclose(sop.first_passage_moment(0), 1.0) 
    assert isclose(sop.probability() * sop.first_passage_moment(1), 1.0)
    # Tests first_passage_probability with eps-downsampling.
    eps=0.9
    sop9 = Syncmer_Parametrized(6,2,ts,eps)
    assert isclose(sop9.first_passage_probability(1), sop.first_passage_probability(1)*(1.0-eps))
    assert isclose(sop.first_passage_moment(0), 1.0) 
    assert isclose(sop.probability() * sop.first_passage_moment(1), 1.0)
    
if __name__ == "__main__":
    main()
