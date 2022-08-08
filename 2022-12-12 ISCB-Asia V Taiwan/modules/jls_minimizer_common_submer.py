#!/usr/bin/env python
"""
Minimizer class.
"""

#from abc import ABC, abstractmethod
import numpy as np
from jls_minimizer import Minimizer # Performs Monte Carlo test.

# S. Schleimer et al. (2003)
# Winnowing: Local Algorithms for Document Fingerprinting
# SIGMOD 2003 76-85

# M.Roberts et al. (2004) 
# Reducing storage requirements for biological sequence comparison
# Bioinformatics 15:3363-3369.

class Minimizer_Common_Submer:
    def __init__(self, w): # (w,k)-minimizer, where k-mer length is irrelevant.
        self.w = w
        self._common = {
            'probability_common_minimizer':dict(), 
            'expectation_common_minimizer':dict(), 
            'common_submer_one_sided':dict(),
        }
    # Returns probability that an alpha-test window contains the same minimizer in two sequences.
    def probability_common_submer(self,alpha): # alpha = test-window length
        #print(alpha)
        w = self.w
        if alpha >= w:
            return 1.0
        # Requires calculation for alpha < w.
        common_submer = dict()
        # Returns the probability_common_minimizer for alpha codes common to two sequences.
        def _common_submer_recursion(a,b,c,d): # w and alpha are passed through calling scope.
            if a+b-2+alpha < w or c+d-2+alpha < w:
                return 0.0
            # Maintains the condition a<=b,c<=d.
            (a,b,c,d) = Minimizer_Common_Submer._canonical_order_4(a,b,c,d) # The symmetry speeds computation.
            # The algorithm assigns the smallest remaining code to an index.
            if (a,b,c,d) not in common_submer:
                ra =  sum([_common_submer_recursion(j,b,c,d) for j in range(1,a)])
                ra += sum([_common_submer_recursion(a,j,c,d) for j in range(1,b)])
                ra += alpha
                ra += sum([_common_submer_recursion(a,b,j,d) for j in range(1,c)])
                ra += sum([_common_submer_recursion(a,b,c,j) for j in range(1,d)])
                common_submer[(a,b,c,d)] = ra/(a+b+c+d-4+alpha)
                #print(a,b,c,d,":",common_submer[(a,b,c,d)])
            return common_submer[(a,b,c,d)]
        # Returns the probability_common_minimizer.
        if alpha not in self._common['probability_common_minimizer']:
            self._common['probability_common_minimizer'][alpha] = _common_submer_recursion(w,w,w,w)
        return self._common['probability_common_minimizer'][alpha]
    # Returns probability that an alpha-test window contains the same minimizer in two sequences.
    def expectation_common_submer(self,alpha): # alpha = test-window length
        w = self.w
        common_submer = dict()
        common_submer_one_sided = self._common['common_submer_one_sided']
        # Returns the expected_common_minimizer_number for alpha codes common to two sequences.
        def _common_submer_recursion(a,b,c,d): # w and alpha are passed through calling scope.
            if a+b-2+alpha < w or c+d-2+alpha < w:
                return 0.0
            (a,b,c,d) = Minimizer_Common_Submer._canonical_order_4(a,b,c,d) # The symmetry speeds computation.
            # The algorithm successively assigns the smallest code to an index.
            if (a,b,c,d) not in common_submer:
                ra =  sum([_common_submer_recursion(j,b,c,d) for j in range(1,a)])
                ra += sum([_common_submer_recursion(a,j,c,d) for j in range(1,b)])
                ra += alpha
                ra += sum([_common_submer_recursion(a,b,j,d) for j in range(1,c)])
                ra += sum([_common_submer_recursion(a,b,c,j) for j in range(1,d)])
                ra += sum([_common_submer_one_sided_recursion(a,c,beta) for beta in range(0,alpha)])
                ra += sum([_common_submer_one_sided_recursion(b,d,beta) for beta in range(0,alpha)])
                common_submer[(a,b,c,d)] = ra/(a+b+c+d-4+alpha)
                #print('common_submer',a,b,c,d,":",common_submer[(a,b,c,d)])
            return common_submer[(a,b,c,d)]
        # Returns the expected_common_minimizer for alpha codes common to two sequences.
        def _common_submer_one_sided_recursion(a,c,beta): # w is passed through calling scope.
            if a-1+beta < w or c-1+beta < w:
                return 0.0
            # Maintains the condition a<=b,c<=d.
            (a,c,beta) = Minimizer_Common_Submer._canonical_order_2(a,c,beta) # The symmetry speeds computation.
            # The algorithm assigns the smallest remaining code to an index.
            if (a,c,beta) not in common_submer_one_sided:
                ra =  sum([_common_submer_one_sided_recursion(j,c,beta) for j in range(1,a)])
                ra += sum([_common_submer_one_sided_recursion(a,j,beta) for j in range(1,c)])
                ra += beta
                ra += sum([_common_submer_one_sided_recursion(a,c,j) for j in range(0,beta)])
                ra += sum([Minimizer.expected_count_of_minimizers(w,j) for j in range(w,beta)])
                common_submer_one_sided[(a,c,beta)] = ra/(a+c-2+beta)
                #print('common_submer_one_sided',a,c,beta,":",common_submer_one_sided[(a,c,beta)])
            return common_submer_one_sided[(a,c,beta)]
        # Returns the probability_common_minimizer.
        if alpha not in self._common['expectation_common_minimizer']:
            self._common['expectation_common_minimizer'][alpha] = _common_submer_recursion(w,w,w,w)
        return self._common['expectation_common_minimizer'][alpha]
    # Returns the canonical order for the arguments (a,b,c,d).
    def _canonical_order_4(a,b,c,d): # (a,b) <= (c,d) and a <= b and c <= d.
        if b < a:
            (a,b,c,d) = (b,a,c,d)
        if d < c:
            (a,b,c,d) = (a,b,d,c)
        if c < a:
            (a,b,c,d) = (c,d,a,b)
        elif c == a and d < b:
            (a,b,c,d) = (c,d,a,b)                
        return (a,b,c,d)
    # Returns the canonical order for the arguments (a,b,c,d).
    def _canonical_order_2(a,c,alpha): # a <= c.
        if c < a:
            (a,c,alpha) = (c,a,alpha)
        return (a,c,alpha)
    # Returns estimates of [probability_common_submer, expectation_common_submer.
    def common_submer_monte_carlo(self, alpha, r): # alpha = test-window length & r = #(realizations) 
        w = self.w
        minimizer = Minimizer(w) # Generates Minimizer realizations.
        css = [0,0] # probability( #(common minimizers) > 1 ) & expectation[ #(common minimizers) ]
        # state consists of indexes [0:w), then [w:w+alpha), then [w+alpha:w+alpha+w).
        for total_realizations in range(r): # Counts the minimizers realized within the sequence.
            common = list(np.random.uniform(size=alpha))
            codes = []
            for i in range(2):
                codes.append(list(np.random.uniform(size=w))+common+list(np.random.uniform(size=w)))
            indicators = [minimizer.indicators(code) for code in codes] # Generates Minimizer realizations.
            is_common_submer = False
            for m in range(w,w+alpha):
                if indicators[0][m] == 1 and indicators[1][m] == 1:
                    if not is_common_submer:
                        css[0] += 1
                        is_common_submer = True
                    css[1] += 1
        return [cs/r for cs in css]

def _test_Minimizer_Common_Submer():
    # Tests probability_common_submer() with a hand calculation.
    minimizer = Minimizer_Common_Submer(2)
    w = minimizer.w
    pcms = [minimizer.probability_common_submer(alpha) for alpha in range(w+1)]
    pcms0 = [0.0, 8/15, 1.0]
    assert np.allclose(pcms,pcms0)
    # Tests probability_common_submer() regressively.
    minimizer = Minimizer_Common_Submer(3)
    w = minimizer.w
    pcms = [minimizer.probability_common_submer(alpha) for alpha in range(w+1)]
    pcms0 = [0.0, 0.3714285714285714, 0.7333333333333334, 1.0]
    assert np.allclose(pcms,pcms0)
    # Tests expectation_common_submer() regressively.
    minimizer = Minimizer_Common_Submer(3)
    w = minimizer.w
    ecms = [minimizer.expectation_common_submer(alpha) for alpha in range(w+1)]
    ecms0 = [0.0, 0.3714285714285714, 0.8066666666666669, 1.2815776815776816]
    assert np.allclose(ecms,ecms0)

def _test_Minimizer_Simulate_Common_Submer():
    r = 10000 #(Monte Carlo realizations)
    np.random.seed(31415)

    minimizer = Minimizer_Common_Submer(3)
    w = minimizer.w
    LENGTH = 2*w
    # Tests common_submer routines against simulation.
    csmc = []
    for alpha in range(LENGTH+1):
        csmc.append(minimizer.common_submer_monte_carlo(alpha, r))
    # Tests probability_common_submer against simulation.
    pcms = [minimizer.probability_common_submer(alpha) for alpha in range(LENGTH+1)]
    pcms_mc = [csmc[alpha][0] for alpha in range(LENGTH+1)]
    assert np.allclose(pcms,pcms_mc,atol=7.0e-03)
    # Tests expectation_common_submer against simulation.
    ecms = [minimizer.expectation_common_submer(alpha) for alpha in range(LENGTH+1)]
    ecms_mc = [csmc[alpha][1] for alpha in range(LENGTH+1)]
    assert np.allclose(ecms,ecms_mc,rtol=2.0e-02)

def main(): 
    _test_Minimizer_Common_Submer()
    _test_Minimizer_Simulate_Common_Submer()
    
if __name__ == "__main__":
    main()
    