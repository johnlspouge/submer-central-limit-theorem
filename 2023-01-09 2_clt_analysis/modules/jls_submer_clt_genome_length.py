#!/usr/bin/env python
"""
Genomic length estimation by Stein's method and a central limit theorem for short-range dependency.
    For input from the abstract class Syncmer jls_syncmer.
"""
# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881
# Returns the syncmer covariances.

from math import sqrt,isinf,isclose
import numpy as np
from scipy.optimize import bisect

from jls_syncmer_parametrized import Syncmer_Parametrized
from jls_submer_clt_util import to_sigma_square_x, to_gamma_x

# Returns the point estimate, the length of the genome.
def length_genome_mean( count_of_submers, mu_y ):
    return count_of_submers/mu_y

# Returns the left endpoint of the domains for f_minus and f_plus.
def lower_limit_L_0( k, sigma_square_y, gamma_y ):
    return max( k, -gamma_y/sigma_square_y)    

# Returns the Wilson score intervals for approximation for L, the genome length.
#    W_Y(L) <= -z_half_alpha
def lambda_minus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ):
    # Returns function f_minus for binary search.
    def _f_minus( L ):
        f_minus = L*mu_y
        f_minus -= z_half_alpha*sqrt(L*sigma_square_y+gamma_y)
        return f_minus - count_of_submers
    L_0 = lower_limit_L_0( k, sigma_square_y, gamma_y )
    L_star = (z_half_alpha/2.0/mu_y)**2*sigma_square_y
    L_star -= gamma_y/sigma_square_y
    L_tilde = max(L_0, L_star)
    if 0.0 <= _f_minus( L_tilde ):
        lambda_minus = [L_0, float('inf')]
    else:
        x = 1.0
        while _f_minus( L_tilde+x ) < 0.0:
            x *= 2.0
        L_minus = bisect( _f_minus, L_tilde, L_tilde+x )
        lambda_minus = [L_minus, float('inf')]
    if L_0 < L_star:
        if 0.0 <= _f_minus( L_0 ):
            L_minus_tilde = bisect( _f_minus, L_0, L_star )
            lambda_minus = [L_0, L_minus_tilde, L_minus, float('inf')]
    assert isinf( lambda_minus[-1] ) 
    return lambda_minus

# Returns the Wilson score intervals for approximation for L, the genome length.
#    W_Y(L) >= z_half_alpha
def lambda_plus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ):
    # Returns function f_plus for binary search.
    def _f_plus( L ):
        f_plus = L*mu_y
        f_plus += z_half_alpha*sqrt(L*sigma_square_y+gamma_y)
        return f_plus - count_of_submers
    L_0 = lower_limit_L_0( k, sigma_square_y, gamma_y )
    if 0.0 < _f_plus( L_0 ):
        lambda_plus = []
    else:
        x = 1.0
        while _f_plus( L_0+x ) <= 0.0:
            x *= 2.0
        L_plus = bisect( _f_plus, L_0, L_0+x )
        lambda_plus = [L_0, L_plus] 
    return lambda_plus

# Returns the Wilson score intervals for approximation for L, the genome length.
#    -z_half_alpha <= W_Y(L) <= z_half_alpha
def lambda_zero( L_0, lambda_minus, lambda_plus ):
    lambda_zero = lambda_minus+lambda_plus
    lambda_zero.sort()
    if lambda_zero[0] == L_0:
        lambda_zero.pop(0)
    else:
        lambda_zero.insert(0, L_0)
    assert isinf( lambda_zero[-1] )
    lambda_zero.pop() # Parity in lambda_zero remains unchanged: it is still interval pairs.
    return lambda_zero

# Returns the standardized_variate_w W(L) in the central limit theorem about the genome length.
def standardized_variate_w(L, submer, count_of_submers): 
    covariances = [ submer.covariance(i) for i in range(submer.get_k()) ]
    return _W( L, count_of_submers, submer.probability(), to_sigma_square_x( covariances ), to_gamma_x( covariances ) )

def _test_standardized_variate_w():
    submer = Syncmer_Parametrized(6, 2, [2])
    covariances = [ submer.covariance(i) for i in range(submer.get_k()) ]
    assert isclose(submer.probability(), 0.2)
    covariances0 = [0.16, -0.04000000000000001, -0.04000000000000001, 0.009999999999999995, 0.004444444444444438, 0.0]
    assert np.allclose(covariances, covariances0)
    assert isclose(to_sigma_square_x( covariances ), 0.028888888888888825)
    assert isclose(to_gamma_x( covariances ), 0.14444444444444457)
    assert isclose(standardized_variate_w(100,submer, 21), 0.5741692517632151)

# Returns the standardized variate from L, the genome length.
def _W( L, count_of_submers, mu_y, sigma_square_y, gamma_y ):
    return (count_of_submers-L*mu_y) / sqrt( L*sigma_square_y+gamma_y )

def _test_lambda():
    def _check(lambda0, value, k, count_of_submers, mu_y, sigma_square_y, gamma_y ):
        for L_est in lambda0:
            if not isinf( L_est ):
                _w = _W( L_est, count_of_submers, mu_y, sigma_square_y, gamma_y )
                assert L_est == lower_limit_L_0( k, sigma_square_y, gamma_y ) or isclose( _w, value ) 
    sigma_square_y = 100.0
    gamma_y = -25.0
    z_half_alpha = 10.0
    mu_y = 5.0
    k = 10
    L_0 = lower_limit_L_0( k, sigma_square_y, gamma_y )
    assert L_0 == 10
    counts_of_submers = [1000000, 10000, 100]
    lambdas = [
        [
            [209146.50210976333, float('inf')],
            [191253.49789023667, 209146.50210976333],
            [10, 191253.49789023667],
        ],[
            [3116.460582894866, float('inf')],
            [1283.539417105132, 3116.460582894866],
            [10, 1283.539417105132],
        ],[
            [438.86068628239263, float('inf')],
            [10, 438.86068628239263],
            [],
        ]
    ]
    for i in range(3):
        count_of_submers = counts_of_submers[i]
        lambda_minus0 = lambda_minus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ) 
        lambda_plus0 = lambda_plus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ) 
        lambda_zero0 = lambda_zero( L_0, lambda_minus0, lambda_plus0 )
        assert np.allclose( lambda_minus0, lambdas[i][0] )
        assert np.allclose( lambda_zero0, lambdas[i][1] )
        assert np.allclose( lambda_plus0, lambdas[i][2] )
        _check(lambda_minus0, -z_half_alpha, k, count_of_submers, mu_y, sigma_square_y, gamma_y)
        _check(lambda_plus0,  +z_half_alpha, k, count_of_submers, mu_y, sigma_square_y, gamma_y)
    k = 200
    L_0 = lower_limit_L_0( k, sigma_square_y, gamma_y )
    assert L_0 == 200
    counts_of_submers = [1000000, 10000, 100]
    lambdas = [
        [
            [209146.50210976333, float('inf')],
            [191253.49789023667, 209146.50210976333],
            [200, 191253.49789023667],
        ],[
            [3116.460582894866, float('inf')],
            [1283.539417105132, 3116.460582894866],
            [200, 1283.539417105132],
        ],[
            [438.86068628239263, float('inf')],
            [200, 438.86068628239263],
            [],
        ]
    ]
    for i in range(3):
        count_of_submers = counts_of_submers[i]
        lambda_minus0 = lambda_minus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ) 
        lambda_plus0 = lambda_plus( count_of_submers, z_half_alpha, k, mu_y, sigma_square_y, gamma_y ) 
        lambda_zero0 = lambda_zero( L_0, lambda_minus0, lambda_plus0 )
        assert np.allclose( lambda_minus0, lambdas[i][0] )
        assert np.allclose( lambda_zero0, lambdas[i][1] )
        assert np.allclose( lambda_plus0, lambdas[i][2] )
        _check(lambda_minus0, -z_half_alpha, k, count_of_submers, mu_y, sigma_square_y, gamma_y)
        _check(lambda_plus0,  +z_half_alpha, k, count_of_submers, mu_y, sigma_square_y, gamma_y)
def _test_lambda_zero():
    assert lambda_zero( 2, [3,4,7,float('inf')], [5,6,9,10] ) == [2, 3, 4, 5, 6, 7, 9, 10]
    assert lambda_zero( 2, [2,3,7,float('inf')], [5,6,9,10] ) == [3, 5, 6, 7, 9, 10]

def main(): 
   _test_standardized_variate_w()
   _test_lambda()
   _test_lambda_zero()
    
if __name__ == "__main__":
    main()
