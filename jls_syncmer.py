#!/usr/bin/env python
"""
Syncmer class
"""

from abc import ABC, abstractmethod
import numpy as np
from jls_submer import Submer

# R. Edgar (2021) 
# Syncmers are more sensitive than minimizers for selecting conserved k‑mers in biological sequences. 
# PeerJ 9: e10805.

class Syncmer(Submer,ABC):
    def __init__(self, k, s):
        super().__init__(k)
        self.s = s
        self.u = k-s
    def get_k(self):
        return super().get_k()
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
    # Returns the syncmer mutation covariances.
    def mutation_covariance(self,i,theta):
        k = self.k
        theta_bar = 1.0-theta
        return self.expectation_product_indicator(i)*(theta_bar**(i+k)-theta_bar**(2*k))
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
    def first_passage_probabilities_monte_carlo(self, r): # r = #(realizations)
        u = self.u
        # Appends to codes until the next syncmer and 
        #   returns the distance to the next syncmer.
        def _distance_to_next_syncmer(codes):
            assert len(codes) == u+1
            i = 0
            while i == 0 or self.indicator(codes[i:i+u+1]) == 0:
                codes.append(np.random.uniform())
                i += 1
            return i
        # The code should work without a window guarantee, so a dictionary is more general than a list.
        fppd = {0:0.0} # first_passage_probabilities conditioned on initial syncmer
        codes = list(np.random.uniform(size=u+1))
        i = _distance_to_next_syncmer(codes)
        codes = codes[i:]
        assert len(codes) == u+1
        assert self.indicator(codes) == 1
        # The codes[0:u] constitute a syncmer.
        for total_realizations in range(r):
            i = _distance_to_next_syncmer(codes)
            codes = codes[i:]
            assert len(codes) == u+1
            assert self.indicator(codes) == 1
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
