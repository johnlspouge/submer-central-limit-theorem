#!/usr/bin/env python
"""
Syncmer_Parametrized class
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
        self.ts = sorted(ts, reverse=True)
        counter = dict(Counter(self.ts))
        for e,v in counter.items():
            assert 0 <= e <= k-s
            assert v == 1 # Every position is unique.
        self.eps = eps
        # Storage of the following probabilities depends on self.eps.
        self._complementary_alpha_run_probability = dict() # _complementary_alpha_run_probability[i]
        # Storage of the following probabilities depends on self.eps.
        self._complementary_alpha_run_probability_mutation = dict() # _complementary_alpha_run_probability[i]
        # Storage of _first_alpha_run_probability does not depend on self.eps.
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
    # Returns the syncmer complementary alpha-run probabilities with downsampling.
    def complementary_alpha_run_probability(self,alpha):
        if alpha <= 0:
            return 1.0
        carp = self.complementary_alpha_run_probability
        _carp = self._complementary_alpha_run_probability
        if alpha not in _carp:
            ts = self.ts
            u = self.u
            # Precomputes for downsampling, as a warm-up for mutation.
            c0 = len(ts)
            min_t = min(ts)
            max_t = max(ts)
            # Computes the downsampling alpha-run probability.
            q = 0.0
            for beta in range(u+alpha):
                if 0 <= beta-max_t and beta-min_t < alpha: # no end-effects
                    c = c0
                else: # end-effects
                    c = len([t for t in ts if 0 <= beta-t < alpha])
                q += self.eps**c*carp(alpha-beta-1)*carp(beta-u)
            q /= u+alpha
            _carp[alpha] = q
        return _carp[alpha]
    # Finds a unique configuration of mutations corresponding to no unmutated syncmers.
    #   Brackets the configuration with (left,right) mutations.
    def first0last0prob_all_mutateds(starts0,k,theta):
        assert 0.0 < theta < 1.0
        # All k-mers overlap.
        starts = sorted(starts0)
        assert starts[0] == 0
        assert starts[-1] < k
        first0last0probs = []
        for last0 in range(k): # mirrored index of last mutation from right starts0[-1]+k-1
            last = starts[-1]+k-1-last0 # index of last mutation
            p0 = theta*(1.0-theta)**last0 
            j0 = len(starts)-1
            while j0 >= 0 and starts[j0] <= last < starts[j0]+k:
                j0 -= 1
            assert j0 == -1 or not starts[j0] <= last < starts[j0]+k # starts[j0] is the last unmutated syncmer.
            if j0 == -1: # All syncmers are mutated.
                for first in range(last+1): # Determines the index of the first mutation.
                    if first < last: # The first & last mutations are distinct.
                        p = p0*theta*(1.0-theta)**first 
                    elif first == last: # The first & last mutations coincide.
                        p = p0*(1.0-theta)**first
                    first0last0probs.append((first,last,p))
            else: # starts[j0] is the last unmutated syncmer, with j0 >= 0.
                for first in range(k):
                    p = p0*theta*(1.0-theta)**first 
                    j = 0
                    while j <= j0 and starts[j] <= first < starts[j]+k:
                        j += 1
                    assert j == j0+1 or not starts[j] <= first < starts[j]+k # starts[j] is the last unmutated syncmer, with j >= 0.
                    if j == j0+1: # No k-mers remain unmutated.
                        first0last0probs.append((first,last,p))
                    elif j <= j0: # starts[j] is unmutated.
                        p *= Syncmer_Parametrized._prob_all_mutated(starts[slice(j,j0+1)],k,theta)
                        first0last0probs.append((first,last,p))
                    else:
                        raise ValueError('j <= j0')
        return sorted(first0last0probs)
    # Returns the total probability that [0:starts[-1]+k) is free of mutation.
    def _prob_all_mutated(starts,k,theta):
        prob = [0.0]*(len(starts)+1)
        prob[0] = 1.0
        i = 0
        while i < len(starts):
            for first in range(k): # first mutation index within starts[i]
                p = theta*(1.0-theta)**first # probability factor from node to subnode
                j = i+1
                while j < len(starts) and starts[j] <= starts[i]+first < starts[j]+k: # The first mutation for starts[i] lies within starts[j].
                    j += 1
                prob[j] += prob[i]*p
            # prob[i] is the total probability of having a mutation in [starts[j]:starts[j]+k) for j in range(i), but not in [starts[i]:starts[-1]+k)
            i += 1
        return prob[-1]
    # Returns the syncmer complementary alpha-run probabilities with downsampling.
    def complementary_alpha_run_probability_mutation(self,alpha,theta):
        if theta == 0.0:
            return Syncmer_Parametrized.complementary_alpha_run_probability(self,alpha)
        if alpha <= 0:
            return 1.0
        print('complementary_alpha_run_probability_mutation',alpha)
        k = self.k
        u = self.u
        ts = self.ts
        # Precomputes (first,last,prob) for full set of syncmers.
        t_min = min(ts)
        t_max = max(ts)
        starts = [t-t_min for t in ts] # 0-offset indexes of syncmers
        flps0 = Syncmer_Parametrized.first0last0prob_all_mutateds(starts,k,theta) # 
        carpm = self.complementary_alpha_run_probability_mutation
        _carpm = self._complementary_alpha_run_probability_mutation
        if alpha not in _carpm:
            ts = self.ts
            u = self.u
            q = 0.0
            for beta in range(u+alpha):
                # Lists the starts of potential syncmers.
                #    For downsampling, not mutation, 
                #        q += self.eps**c*carpm(beta-u)*carpm(alpha-beta-1)
                if 0 <= beta-t_max and beta-t_min < alpha: # SNAFU: All syncmers are present.
                    for flp in flps0:
                        first,last,prob = flp
                        lt = beta-u+min(first-t_max,0)
                        rt = alpha-beta-1-max(last-t_max,0)
                        q += prob*carpm(lt,theta)*carpm(rt,theta)
                else:
                    ts0 = [t for t in ts if 0 <= beta-t < alpha]
                    if not ts0:
                        q += carpm(alpha-beta-1,theta)*carpm(beta-u,theta)
                    else:
                        t_min0 = min(ts0)
                        t_max0 = max(ts0)
                        starts0 = [t-t_min0 for t in ts0]
                        flps00 = Syncmer_Parametrized.first0last0prob_all_mutateds(starts0,k,theta)
                        for flp in flps00:
                            first,last,prob = flp
                            lt = beta-u+min(first-t_max0,0)
                            rt = alpha-beta-1-max(last-t_max0,0)
                            q += prob*carpm(lt,theta)*carpm(rt,theta)
            q /= u+alpha
            _carpm[alpha] = q
        return _carpm[alpha]
    # Returns the syncmer complementary alpha-run probabilities.
    def alpha_run_probability(self,alpha):
        return 1.0-self.complementary_alpha_run_probability(alpha)
    # Returns the syncmer complementary alpha-run probabilities.
    def first_passage_probability(self,alpha,theta=0.0):
        if alpha == 0:
            return 0.0
        carp = self.complementary_alpha_run_probability_mutation
        fpp = 0.0
        fpp += carp(alpha+1,theta)
        fpp += -2.0*carp(alpha,theta)
        fpp += carp(alpha-1,theta)
        fpp /= self.probability()*(1.0-theta)**self.k
        return fpp
# Tests parametrized syncmers (closed syncmers).
def _test_Syncmer_Closed():
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
# Tests parametrized syncmers (open syncmers).
def _test_Syncmer_Open():
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

def _test_prob_all_mutateds():
    first0last0prob_all_mutateds0 = [(0, 3, 0.00081), (0, 4, 0.0081), (0, 5, 0.0809919), (0, 6, 0.8098461), (1, 3, 8.1e-05), (1, 4, 0.00081), (1, 5, 0.0081), (1, 6, 0.0809919), (2, 3, 8.1e-06), (2, 4, 8.1e-05), (2, 5, 0.00081), (2, 6, 0.0081), (3, 3, 9.0e-07), (3, 4, 8.1e-06), (3, 5, 8.1e-05), (3, 6, 0.00081)]
    assert isclose(Syncmer_Parametrized._prob_all_mutated([0,1,2,3],4,0.9),0.99963)
    assert isclose(sum([c[2] for c in first0last0prob_all_mutateds0]),0.99963)
    assert np.allclose(Syncmer_Parametrized.first0last0prob_all_mutateds([0,1,2,3],4,0.9), first0last0prob_all_mutateds0)
    first0last0prob_all_mutateds0 = [(0, 3, 0.00081), (0, 4, 0.0081), (0, 5, 0.081), (0, 6, 0.809919), (1, 3, 8.1e-05), (1, 4, 0.00081), (1, 5, 0.0081), (1, 6, 0.0809919), (2, 3, 8.1e-06), (2, 4, 8.1e-05), (2, 5, 0.00081), (2, 6, 0.0081), (3, 3, 9.0e-07), (3, 4, 8.1e-06), (3, 5, 8.1e-05), (3, 6, 0.00081)]
    assert np.allclose(Syncmer_Parametrized.first0last0prob_all_mutateds([0,2,3],4,0.9), first0last0prob_all_mutateds0)
    assert isclose(Syncmer_Parametrized._prob_all_mutated([0,2,3],4,0.9), 0.999711)
    assert isclose(sum([c[2] for c in first0last0prob_all_mutateds0]), 0.999711)
    first0last0prob_all_mutateds0 = [(0, 0, 0.09), (0, 1, 0.81), (1, 1, 0.09)]
    assert np.allclose(Syncmer_Parametrized.first0last0prob_all_mutateds([0],2,0.9), first0last0prob_all_mutateds0)
    assert isclose(Syncmer_Parametrized._prob_all_mutated([0],2,0.9),0.99)
    assert isclose(sum([c[2] for c in first0last0prob_all_mutateds0]),0.99)
    first0last0prob_all_mutateds0 = [(0, 1, 0.081), (0, 2, 0.81), (1, 1, 0.009), (1, 2, 0.081)]
    assert np.allclose(Syncmer_Parametrized.first0last0prob_all_mutateds([0,1],2,0.9), first0last0prob_all_mutateds0)
    assert isclose(Syncmer_Parametrized._prob_all_mutated([0,1],2,0.9),0.981)
    assert isclose(sum([c[2] for c in first0last0prob_all_mutateds0]),0.981)
    first0last0prob_all_mutateds0 = [(0, 2, 0.0081), (0, 3, 0.081), (0, 4, 0.81), (1, 2, 0.00081), (1, 3, 0.0081), (1, 4, 0.081), (2, 2, 9.0e-05), (2, 3, 0.00081), (2, 4, 0.0081)]
    assert np.allclose(Syncmer_Parametrized.first0last0prob_all_mutateds([0,2],3,0.9), first0last0prob_all_mutateds0)
    assert isclose(Syncmer_Parametrized._prob_all_mutated([0,2],3,0.9),0.99801)
    assert isclose(sum([c[2] for c in first0last0prob_all_mutateds0]),0.99801)

def _test_Syncmer_Parametrized_Simulate():
    r = 100000 #(Monte Carlo realizations)
    np.random.seed(31415)
    k = 6
    s = 2
    ts = [2,4]
    sp = Syncmer_Parametrized(k,s,ts)
    """
    # Tests first_passage_probabilities by simulation.
    fpps_mc = sp.first_passage_probabilities_monte_carlo(r)
    assert isclose(sum(fpps_mc),1.0)
    fpps_mc0 = [0.0, 0.16633, 0.47696, 0.16963, 0.11318, 0.04285, 0.02005, 0.00728, 0.00248, 0.00092, 0.00024, 5e-05, 3e-05]
    assert np.allclose(fpps_mc, fpps_mc0)
    # The first_passage_probability inverts Shaw's formula to verify the alpha_run_probability.
    fpps = [ sp.first_passage_probability(alpha) for alpha,f in enumerate(fpps_mc) ]
    fpps0 = [0.0, 0.16666666666666663, 0.47619047619047616, 0.1696428571428571, 0.11342592592592592, 0.04312169312169312, 0.019859307359307354, 0.00709726230559564, 0.002688322479989146, 0.0008857544571830286, 0.0002945030425189155, 8.967000720472941e-05, 2.7017391228700758e-05]
    assert np.allclose(fpps, fpps0)
    # Tests first_passage_probabilities with downsampling eps by simulation (no mutation theta=0.0).
    eps = 0.1
    theta = 0.0
    sp = Syncmer_Parametrized(k,s,ts,eps)
    fpps_mc = sp.first_passage_probabilities_monte_carlo(r, eps, theta)
    assert isclose(sum(fpps_mc),1.0)
    fpps_mc0 = [0.0, 0.15064, 0.43168, 0.16972, 0.11761, 0.06271, 0.03521, 0.01672, 0.00811, 0.00419, 0.00188, 0.00076, 0.00037, 0.00021, 0.00013, 6e-05]
    assert np.allclose(fpps_mc,fpps_mc0)
    # The first_passage_probability inverts Shaw's formula to verify the alpha_run_probability with downsampling.
    fpps = [ sp.first_passage_probability(alpha) for alpha,f in enumerate(fpps_mc) ]
    fpps0 = [0.0, 0.1500000000000001, 0.43071428571428544, 0.1709196428571429, 0.11833660714285725, 0.06221193749999994, 0.03518347134740263, 0.01682962781858765, 0.008410472303618256, 0.003929942195703741, 0.0018637777932891712, 0.0008581962968659726, 0.0003990719989349006, 0.00018384035198386816, 8.527432968391637e-05, 3.9488892042046585e-05]
    assert np.allclose(fpps,fpps0)
    """
    # Tests first_passage_probabilities with downsampling eps by simulation (no mutation theta=0.0).
    eps = 0.0
    theta = 0.1
    sp = Syncmer_Parametrized(k,s,ts,eps)
    fpps_mc = sp.first_passage_probabilities_monte_carlo(r, eps, theta)
    assert isclose(sum(fpps_mc),1.0)
    fpps_mc0 = [
        0.0, 0.14662, 0.38534, 0.12376, 0.08096, 
        0.05535, 0.05091, 0.03921, 0.02865, 0.02184, 
        0.01649, 0.01231, 0.00968, 0.00717, 0.00538, 
        0.00393, 0.00295, 0.00237, 0.00174, 0.00143, 
        0.00099, 0.0007, 0.00061, 0.00033, 0.00034, 
        0.00029, 0.00017, 0.00015, 8e-05, 7e-05, 
        5e-05, 3e-05, 3e-05, 0.0, 0.0, 
        4e-05, 1e-05, 0.0, 2e-05
    ]
    fpps_mc0 = [sp.first_passage_probability(alpha,theta) for alpha,f in enumerate(fpps_mc)]
    fpps_mc00 = [
        0.0, 0.14589016666666638, 0.3856767650885492, 0.13365188294255867, 0.07823208656871342, 
        0.020205643591354978, 0.013210708225969802, 0.026108351296338955, 0.03122359349373644, 0.026234695652148806, 
        0.020835082055570903, 0.016876047818387897, 0.014307249374428642, 0.01252709711292685, 0.010871000572336663, 
        0.009288655694309212, 0.007905878219766784, 0.006747003426805521, 0.005778266810175532, 0.004955279193706475, 
        0.004245957714682454, 0.003634141801785798, 0.0031094498126601317, 0.002661066157869699, 0.0022779566821390833, 
        0.0019502043783575917, 0.0016695357946798315, 0.0014291635809622228, 0.0012233697716326271, 0.0010472206153486808, 
        0.0008964485458523153, 0.0007673877437382648, 0.0006569056862584064, 0.0005623276545504647, 0.00048136593766193473, 
        0.0004120611517604218, 0.0003527349200168107, 0.0003019502371449777, 0.00025847716818792253]
    """
    #fpps_mc0 = [0.0, 0.15064, 0.43168, 0.16972, 0.11761, 0.06271, 0.03521, 0.01672, 0.00811, 0.00419, 0.00188, 0.00076, 0.00037, 0.00021, 0.00013, 6e-05]
    #assert np.allclose(fpps_mc,fpps_mc0)
    """

def main():
    _test_prob_all_mutateds()
    _test_Syncmer_Closed()
    _test_Syncmer_Open()
    """
    _test_Syncmer_Parametrized_Simulate()
    """
    
if __name__ == "__main__":
    main()
