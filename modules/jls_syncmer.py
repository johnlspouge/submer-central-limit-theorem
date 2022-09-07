#!/usr/bin/env python
"""
Syncmer class
"""

from abc import ABC, abstractmethod
import numpy as np
from scipy.stats import bernoulli
from collections import deque
from jls_submer import Submer

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.

class Syncmer(Submer,ABC):
    def __init__(self, k, s):
        super().__init__(k)
        self.s = s
        self.u = k-s
    @abstractmethod
    # Returns the indicator corresponding to a list of (u+1) s-codes.
    def indicator(self, codes):
        pass
    # Returns a list of indicators for the indexes of minimizers in a list of codes.
    #   len(indicators) == len(codes).
    def indicators(self, codes):
        u = self.u
        indicators = [] 
        for i in range(len(codes)-u): # c-(k-s) k-mers yield c s-codes.
            indicators.append(self.indicator(codes[i:i+u+1]))
        return indicators
    # Returns a complete set of indicators.
    #   Updates the indicators by writing 1 at all syncmers after last.
    #   Assumes the last syncmer computed is codes[last].
    #   Assumes len(indicators) == len(codes)-(k-s).
    def update_indicators(self,indicators,codes,last):
        u = self.u
        assert len(indicators) == len(codes)-u 
        assert last < len(indicators)
        assert indicators[last] == 1
        last_new = last
        for i in range(last+1,len(codes)-u):
            indicator = self.indicator(codes[i:i+u+1])
            indicators[i] = indicator
            if indicator == 1:
                last_new = i
        return last_new
    # Returns estimates of expectation_product_indicators as a list of length u+1,
    #     because E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > u.
    def expectation_product_indicators_monte_carlo(self, r): # r = #(realizations) 
        u = self.u
        epis = [0]*(u+1) # expectation_product_indicators
        codes = np.random.uniform(size=2*(u+1))
        indicators = self.indicators(codes)
        for i in range(r):
            i0 = indicators[0]
            epis = [a+b*i0 for (a,b) in zip(epis,indicators)] # expectation_product_indicators
            codes = np.random.uniform(size=2*(u+1))
            indicators = self.indicators(codes)
        epis = [a/r for a in epis]
        return epis
    # Returns a list indexed range(n) of non-0 estimates of first_passage_probabilities[0...n).
    #   tested in jls_syncmer_parametrized for downsampling with 0.0 < eps.
    #   untested for mutation with 0.0 < theta.
    def first_passage_probabilities_monte_carlo(self, 
            r, # #(realizations)
            eps=0.0, # downsampling rejection prob
            theta=0.0): # mutation prob
        assert 0.0 <= eps < 1.0
        assert 0.0 <= theta < 1.0
        k = self.k
        u = self.u
        # Appends to codes until the next syncmer and returns its distance.
        def _distance_to_next_syncmer(codes, eps, theta, q):
            assert len(codes) == u+1
            assert theta == 0.0 or len(q) == k
            i = 0
            while True:
                codes.append(np.random.uniform())
                i += 1
                if self.indicator(codes[i:i+u+1]) == 0: # The k-mer is not a syncmer.
                    continue 
                if 0.0 < eps and bernoulli.rvs(eps): # The syncmer is downsampled.
                    continue
                if 0.0 < theta: # The syncmer can be mutated.
                    n = min(i,k)
                    bs = bernoulli.rvs(theta,size=n)
                    for j in range(n):
                        q.pop()
                        q.appendleft(bs[j])
                    assert len(q) == k
                    if not all(b == 0 for b in q): # The syncmer is mutated.
                        continue
                break
            return i
        q = [0]*k # mutation status
        if 0.0 < theta: 
            q = deque(bernoulli.rvs(theta, size=k))
        # The code needs to work without a window guarantee, so a dictionary is more general than a list.
        fppd = {0:0} # first_passage_probabilities conditioned on initial syncmer
        codes = list(np.random.uniform(size=u+1))
        # Initializes codes[0:u] at the next syncmer.
        i = _distance_to_next_syncmer(codes, eps, theta, q)
        codes = codes[i:]
        assert len(codes) == u+1
        assert self.indicator(codes) == 1
        if 0.0 < theta: # mutation
            assert all(b == 0 for b in q)
        # The codes[0:u] are an unmutated syncmer.
        for total_realizations in range(r):
            i = _distance_to_next_syncmer(codes, eps, theta, q)
            codes = codes[i:]
            assert len(codes) == u+1
            assert self.indicator(codes) == 1
            if 0.0 < theta: # mutation
                assert all(b == 0 for b in q)
            if i not in fppd:
                fppd[i] = 1
            else:
                fppd[i] += 1
        # Converts a dictionary to a list by padding internal 0s.
        def _to_list(dictionary):
            fpps = []
            for k in range(max(dictionary.keys())+1):
                if k in dictionary:
                    fpps.append(dictionary[k])
                else:
                    fpps.append(0)
            return fpps
        fpps = _to_list(fppd)
        assert sum(fpps) == r
        return [fpp/r for fpp in fpps]

def main(): 
    pass
    
if __name__ == "__main__":
    main()
