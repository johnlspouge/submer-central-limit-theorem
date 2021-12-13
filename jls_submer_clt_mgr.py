#!/usr/bin/env python
"""
Manager routines for Stein's method and a central limit theorem for short-range dependency in submers producing genome length.
"""

# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881

from math import isclose
from numpy import allclose

from jls_syncmer_closed import Syncmer_Closed
from jls_syncmer_open import Syncmer_Open
from jls_minimizer import Minimizer
from jls_syncmer_clt_util import to_sigma_square_x, to_gamma_x, to_delta, z_half_alpha 
from jls_submer_clt_genome_length import length_genome_mean, lower_limit_L_0, lambda_minus, lambda_plus, lambda_zero
from jls_syncmer_clt_theta import theta_all, theta_mean

# Returns the Wilson score intervals for genomic length L:
#   "left": Pr(L in the left interval) <= alpha_0 / 2
#   "right": Pr(L in the right interval) <= alpha_0 / 2
#   "confidence": Pr(L in the confidence interval) >= 1.0 - alpha_0

# The following 3 subroutines for submers (closed syncmers, open syncmers, or minimizers) have the following arguments:
#   count_of_submers: the (possibly unassembled) counts for the submers
#   alpha0: the upper bound on the size of the test for left, right, or confidence intervals
#   (k,s), (k,s,t), or (w,k) are the parameters governing closed syncmers, open syncmers, or minimizers, where k is the k-mer size of the submer. 
#   is_estimate: (None => qualitative method) or (True => estimate method)

# Returns None for the estimate method if the bound in Stein's method is too large to be useful.

# Syncmer_Closed(k, s) # k-mer syncmers with s-codes
def to_syncmer_closed_Wilson_score_intervals_for_length( count_of_submers, alpha_0, k, s, is_estimate=None ):
    syncmer_closed = Syncmer_Closed(k, s)
    return _to_submer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, syncmer_closed, is_estimate )   

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes
def to_syncmer_open_Wilson_score_intervals_for_length( count_of_submers, alpha_0, k, s, t, is_estimate=None ):
    syncmer_open = Syncmer_Open(k, s, t)
    return _to_submer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, syncmer_open, is_estimate )
    
# Minimizer(w, k) # (w,k)-minimizers = w-windows of k-mers
def to_minimizer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, w, k, is_estimate=None ):
    minimizer = Minimizer(w, k)
    return _to_submer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, minimizer, is_estimate )
    
# common Wilson score subroutine for submers
#   Returns None if Stein's method fails to produce a useful bound.
def _to_submer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, submer, is_estimate=None ):
    k = submer.get_k()
    mu_y = submer.probability()
    covariances = [ submer.covariance(i) for i in range(k) ]
    if is_estimate is None:
        delta = 0.0 # qualitative method
    else:
        L_hat = length_genome_mean( count_of_submers, mu_y ) 
        delta = to_delta( covariances, L_hat, mu_y, mu_y ) # estimate method
        if alpha_0 <= delta:
            return None # Returns None for the estimate method if the bound in Stein's method is too large to be useful.
    z_half_alpha_ = z_half_alpha( alpha_0, delta )
    sigma_square_y = to_sigma_square_x( covariances )
    gamma_y = to_gamma_x( covariances )
    L_0 = lower_limit_L_0( k, sigma_square_y, gamma_y )
    left = lambda_minus( count_of_submers, z_half_alpha_, k, mu_y, sigma_square_y, gamma_y )
    right = lambda_plus( count_of_submers, z_half_alpha_, k, mu_y, sigma_square_y, gamma_y )
    confidence = lambda_zero( L_0, left, right )
    return {
        "left":left,
        "confidence":confidence,
        "right":right,
    }   

# Returns the point estimate, the length of the genome.
def length_point_estimate( count_of_submers, mu_y ): # mu_y is the submer density, submer.probability()
    return length_genome_mean( count_of_submers, mu_y )

def _test_to_submer_Wilson_score_intervals_for_length():
    k = 10
    s = 3
    syncmer_closed = Syncmer_Closed(k,s)
    counts_of_submers = [100, 10000, 1000000]
    lambdas0 = [
       {'left': [435.4890933863735, float('inf')], 'confidence': [367.3508315817253, 435.4890933863735], 'right': [10, 367.3508315817253]},
       {'left': [40338.49819437848, float('inf')], 'confidence': [39664.3417305896, 40338.49819437848], 'right': [10, 39664.3417305896]},
       {'left': [4003371.839017747, float('inf')], 'confidence': [3996631.000907222, 4003371.839017747], 'right': [10, 3996631.000907222]},
    ]
    for i in range(3):
        count_of_submers = counts_of_submers[i]
        assert isclose(2/(k-s+1), 0.25) # syncmer density
        assert isclose(length_point_estimate( count_of_submers, syncmer_closed.probability() ), count_of_submers / 0.25)
        alpha_0 = 0.05
        is_estimate=None # In the estimate method, the bound in Stein's method was always too large to be useful at alpha_0 = 0.05.
        d = to_syncmer_closed_Wilson_score_intervals_for_length( count_of_submers, alpha_0, k, s, is_estimate ) # is_estimate=None is the default.
        for key in d:
            assert allclose( d[key], lambdas0[i][key] )

# Returns the Wilson score intervals for theta, mutation probability per base:
#   "left": Pr(L in the left interval) <= alpha_0 / 2
#   "right": Pr(L in the right interval) <= alpha_0 / 2
#   "confidence": Pr(L in the confidence interval) >= 1.0 - alpha_0

# The following 2 subroutines for submers (closed syncmers or open syncmers) have the following arguments:
#   count_of_submers: the (possibly unassembled) counts for the submers
#   count_of_unmutated_syncmers: the (possibly unassembled) counts for the submers common to the 2 sequences under scrutiny
#   alpha0: the upper bound on the size of the test for left, right, or confidence intervals
#   (k,s) or (k,s,t) are the parameters governing closed syncmers or open syncmers, where k is the k-mer size of the submer. 
#   length_genome: (None => estimate method) or (integer => actual genome length from data)

# Syncmer_Closed(k, s) # k-mer syncmers with s-codes
def to_syncmer_closed_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, k, s, length_genome=None ):
    syncmer_closed = Syncmer_Closed(k, s)
    return _to_syncmer_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, syncmer_closed, length_genome )

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes
def to_syncmer_open_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, k, s, t, length_genome=None ):
    syncmer_open = Syncmer_Open(k, s, t)
    return _to_syncmer_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, syncmer_open, length_genome )
    
def _to_syncmer_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, syncmer, length_genome=None ):
    if length_genome is None:
        mu_y = syncmer.probability()
        L_hat = length_point_estimate( count_of_syncmers, mu_y ) # estimated genome length
    elif not isinstance(length_genome, int) or length_genome <= 0:
        raise Exception("'length_genome' must be a positive integer.")
    else:
        L_hat = length_genome # actual genome length
    delta = 0.0 # The standard deviation in the denominator of the standardized variate can be arbitrarily small, frustrating Stein's method.
    z_half_alpha_ = z_half_alpha( alpha_0, delta )
    thetas = theta_all( count_of_syncmers, count_of_unmutated_syncmers, z_half_alpha_, syncmer, L_hat )
    return {
        "left":thetas[0],
        "confidence":thetas[1],
        "right":thetas[2],
    }   

# Returns the point estimate, the mutation probability per letter.
def theta_point_estimate( count_of_syncmers, count_of_unmutated_syncmers, k ):
    return theta_mean( count_of_syncmers, count_of_unmutated_syncmers, k )

def _test_to_syncmer_Wilson_score_intervals_for_theta():
    k = 10
    s = 3
    #syncmer_closed = Syncmer_Closed(k,s)
    counts_of_syncmers = [100, 10000, 1000000]
    thetas0 = [
        {'left': [0.0, 0.0038246028859374994], 'confidence': [0.0038246028859374994, 0.026768396719205736], 'right': [0.026768396719205736, 1.0]},
        {'left': [0.0, 0.009450735795312502], 'confidence': [0.009450735795312502, 0.011614293370478451], 'right': [0.011614293370478451, 1.0]},
        {'left': [0.0, 0.010373162301562497], 'confidence': [0.010373162301562497, 0.01058935753615523], 'right': [0.01058935753615523, 1.0]},
    ]
    for i in range(3):
        count_of_syncmers = counts_of_syncmers[i]
        count_of_unmutated_syncmers = count_of_syncmers * 0.9
        assert isclose(theta_point_estimate( count_of_syncmers, count_of_unmutated_syncmers, k ), 0.010480741793785553) # The point estimate depends only on syncmer ratios.
        alpha_0 = 0.05
        d = to_syncmer_closed_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, k, s, length_genome=None ) 
        for key in d:
            assert allclose( d[key], thetas0[i][key] )
        
def main(): 
    _test_to_submer_Wilson_score_intervals_for_length()
    _test_to_syncmer_Wilson_score_intervals_for_theta()
    
if __name__ == "__main__":
    main()
