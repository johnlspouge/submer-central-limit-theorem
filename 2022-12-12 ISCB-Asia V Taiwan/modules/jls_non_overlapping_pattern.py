#!/usr/bin/env python
"""
non-overlapping pattern class
"""
# M.C. Frith, L. Noé, G. Kucherov (2020) 
# Minimally-overlapping words for sequence similarity search. 
# Bioinformatics 36: 5344-50.

from math import isclose
import numpy as np
import random
import re
from jls_pattern_util import is_dna_iupac, to_regex, to_probability, overlap
from jls_non_overlapping_pattern_prob import Non_Overlapping_Pattern_Prob

class Non_Overlapping_Pattern(Non_Overlapping_Pattern_Prob): 
    # Argument patterns is a non-overlapping pattern or a list of non-overlapping patterns
    #   Patterns are strings of (possibly ambiguous) IUPAC DNA nucleotides.
    #   Default search is case-insensitive (regex_flags=re.IGNORECASE).
    def __init__(self, 
                 k, # k-mer length 
                 patterns, # a IUPAC string or list of IUPAC strings
                 nucleotide2probability = {'a':0.25,'t':0.25,'c':0.25,'g':0.25}, # independent letters model : a dict keys=chars a|c|g|t values=probabilities, where sum(probabilities) == 1.0
                 regex_flags=re.IGNORECASE): # flags for sequence matching (re.IGNORECASE)
        if isinstance(patterns,str):
            patterns = [patterns]
        self._patterns = patterns
        self._nucleotide2probability = nucleotide2probability
        self._check()
        regex = '^' # The patterns are submer prefixes.
        regex += to_regex(patterns)
        self._regex = re.compile(regex, flags=regex_flags)
        length2probability = self.to_length2probability()
        super().__init__(k, length2probability)
    def _check(self):
        patterns = self._patterns
        nucleotide2probability = self._nucleotide2probability
        # Checks patterns.
        for pattern in patterns:
            if not is_dna_iupac(pattern):
                raise Exception(f'The pattern "{pattern}" is not a IUPAC DNA string.')
        overlapc = overlap(patterns)
        if overlapc:
            raise Exception(f'The patterns "{patterns}" have non-empty overlap "{overlapc}".')
        # Checks nucleotides.
        nucleotides = ''.join(sorted(nucleotide2probability.keys())).lower()
        if nucleotides != 'acgt':
            raise Exception('The list of nucleotide2probability.keys() is not the 4 unambiguous DNA nucleotides.')
        probs = nucleotide2probability.values()
        if not all([ 0.0 <= prob for prob in probs ]):
            raise Exception('The list of nucleotide2probability.values() has negative elements: {probs}.')
        if not isclose(sum(probs), 1.0):
            raise Exception('The list of nucleotide2probability.values() do not sum to 1.0: {probs}.')
    # In an independent letters model, nucleotide2probability gives the probability of each of the 4 unambiguous nucleotide.
    def to_length2probability(self):
        length2probability = {}
        for pattern in self._patterns:
            length = len(pattern)
            prob = np.prod([to_probability(iupac, self._nucleotide2probability) for iupac in pattern])
            length2probability[length] = length2probability.get(length, 0.0) + prob
        return length2probability
    # Returns pattern indicator for k-mer.
    #   Assumes len(kmer_letters == self.k).
    def indicator(self, kmer_letters):
        return bool(re.search(self._regex, kmer_letters))
    # Returns a list of pattern indicators for string.
    #   len(indicators) == len(string - k + 1).
    def indicators(self, string):
        k = self.k
        indicators = []
        i = 0
        while i < len(string)-k+1:
            indicators.append(self.indicator(string[i:i+k]))
            i += 1
        return indicators
    # Returns estimates of expectation_product_indicators as a list of length u+1,
    #     because E[ Y[0]*Y[i] ] = E[ Y[0] ]*E[ Y[i] ] for i > u.
    def expectation_product_indicators_monte_carlo(self, r): # r = #(realizations) 
        k = self.k
        nucleotides = list(self._nucleotide2probability.keys())
        probabilities = list(self._nucleotide2probability.values())
        epis = [0]*(k+1) # expectation_product_indicators
        letters = random.choices(nucleotides, weights=probabilities, k=2*self.k)
        letters = ''.join(letters)
        indicators = self.indicators(letters)
        assert len(indicators) == self.k + 1
        for i in range(r):
            i0 = indicators[0]
            epis = [a+b*i0 for (a,b) in zip(epis,indicators)] # expectation_product_indicators
            letters = random.choices(nucleotides, weights=probabilities, k=2*self.k)
            letters = ''.join(letters)
            indicators = self.indicators(letters)
        epis = [a/r for a in epis]
        return epis
    # Returns a list indexed range(n) of non-0 estimates of first_passage_probabilities[0...n).
    def first_passage_probabilities_monte_carlo(self, r): # r = #(realizations)
        k = self.k
        nucleotides = list(self._nucleotide2probability.keys())
        probabilities = list(self._nucleotide2probability.values())
        # Appends to letters until the next submer and 
        #   returns the letters from present k-mer to the next submer.
        def _to_next_submer(letters):
            assert len(letters) == k
            i = 0
            while i == 0 or self.indicator(letters[i:i+k]) == 0:
                letters += random.choices(nucleotides, weights=probabilities, k=1)[0]
                i += 1
            return letters
        # The code should work without a window guarantee, so a dictionary is more general than a list.
        fppd = {0:0.0} # first_passage_probabilities conditioned on initial syncmer
        letters = random.choices(nucleotides, weights=probabilities, k=self.k)
        letters = ''.join(letters)
        letters = _to_next_submer(letters)
        i = len(letters) - k # letters added
        letters = letters[i:]
        assert len(letters) == k
        assert self.indicator(letters) == 1
        # The codes[0:k) constitute a syncmer.
        for total_realizations in range(r):
            letters = _to_next_submer(letters)
            i = len(letters) - k # letters added
            letters = letters[i:]
            assert len(letters) == k
            assert self.indicator(letters) == 1
            fppd[i] = fppd.get(i, 0) + 1
        # Converts a dictionary to a list by padding internal 0s.
        def _to_list(dictionary):
            fpps = []
            for k in range(max(dictionary.keys())+1):
                fpps.append(dictionary.get(k, 0))
            return fpps
        fpps = _to_list(fppd)
        assert sum(fpps) == r
        return [fpp/r for fpp in fpps]

def _test_Non_Overlapping_Pattern():
    patterns = ['a', 'cg']
    # Tests indicator() and indicators() by hand calculation.
    nop = Non_Overlapping_Pattern(3,patterns)
    assert nop.indicator('aaa') == 1
    assert nop.indicator('att') == 1
    assert nop.indicator('taa') == 0
    assert nop.indicator('cga') == 1
    assert nop.indicator('tcg') == 0
    assert nop.indicators('aaaa') == [True, True]
    assert nop.indicators('atta') == [True, False]
    assert nop.indicators('taaa') == [False, True]
    assert nop.indicators('cgaa') == [True, False]
    assert nop.indicators('tcga') == [False, True]
    # Tests lagging moments for island calculations.
    k=6
    patterns = ['a', 'cg']
    nucleotide2probability = {'a':0.15,'c':0.2,'g':0.5,'t':0.15}
    nop = Non_Overlapping_Pattern(k,patterns,nucleotide2probability)
    p = nop.probability() # syncmer density
    assert isclose(nop.first_passage_moment(0,k)*p, 0.03962409765625) 
    assert isclose(nop.first_passage_moment(1,k)*p, 0.136140859375)
    assert isclose(nop.first_passage_moment(2,k)*p, 0.799362265625) 
    # f=[0.0, 0.15, 0.2275, 0.178375, 0.12886875, 0.0917009375, 0.065058921875, 0.04612998984375, 0.032704599179687495, 0.02318591031835937, 0.016437563852636717, 0.011653338242905273, 0.008261581121205811, 0.005857010128734412, 0.004152300497303669, 0.0029437544098346773, 0.002086961198629109, 0.0014795415778512749, 0.0010489142213106727, 0.0007436229303289443, 0.0005271880686485354, 0.0003737475653183606, 0.00026496662365575296, 0.00018784687357555394]    
def _test_Non_Overlapping_Pattern_Simulate():
    r = 100000 #(Monte Carlo realizations) # tested for 100*r, too
    np.random.seed(31415)
    non_overlapping_pattern = Non_Overlapping_Pattern(4,'ryy',
                                {'a':0.1,'t':0.4,'c':0.2,'g':0.3})
    nop = non_overlapping_pattern
    k = nop.k
    # Tests expectation_product_indicators by simulation.
    epis_mc = nop.expectation_product_indicators_monte_carlo(r)
    epis = [ nop.expectation_product_indicator(i) for i in range(len(epis_mc)) ] # length k
    assert np.allclose( epis_mc, epis, atol=2.0e-03 )
    assert len(epis_mc) == k+1 # expectation_product_indicators_monte_carlo return len(list) == k.
    # Tests first_passage_probabilities by simulation.
    fpps_mc = nop.first_passage_probabilities_monte_carlo(r)
    fpps = [ nop.first_passage_probability(i) for i in range(len(fpps_mc)) ]
    assert np.allclose( fpps_mc, fpps, atol=3.0e-03 )
    assert isclose(sum(fpps_mc),1.0)

def main(): 
    _test_Non_Overlapping_Pattern()
    _test_Non_Overlapping_Pattern_Simulate()
    
if __name__ == "__main__":
    main()
