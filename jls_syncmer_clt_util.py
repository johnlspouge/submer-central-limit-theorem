#!/usr/bin/env python
"""
Utility routines for Stein's method and a central limit theorem for short-range dependency.
    For input from the abstract class Syncmer jls_syncmer.
"""
# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881

from math import sqrt,pi,isclose
from scipy.stats import norm

# Returns the syncmer sigma_square_x.
def to_sigma_square_x( covariances ):
    k = len(covariances)
    sigma_square_x = 0.0
    for i in range(1,k):
        sigma_square_x += covariances[i]
    sigma_square_x *= 2.0
    sigma_square_x += covariances[0]
    return sigma_square_x
def _test_to_sigma_square_x():
    assert to_sigma_square_x( [2,3,4,5] ) == 26

# Returns the syncmer gamma_x.
def to_gamma_x( covariances ):
    k = len(covariances)
    gamma_x = 0.0
    for i in range(1,k):
        gamma_x += i*covariances[i]
    gamma_x *= -2.0
    return gamma_x
def _test_to_gamma_x():
    assert to_gamma_x( [2,3,4,5] ) == -52

# Returns the standard deviation sigma_l for central limit theorems with short-range dependency.
def to_sigma_l( covariances, length_genome ):
    if length_genome < len(covariances):
        raise Exception("length_genome < len(covariances)")
    variance_l = length_genome*to_sigma_square_x( covariances )+to_gamma_x( covariances )
    try:
        sigma_l = sqrt(variance_l)
    except:
        raise Exception("sqrt(variance_l) : variance_l < 0.0")
    return sigma_l
def _test_to_sigma_l():
    assert isclose( to_sigma_l( [2,3,4,5], 100 ), sqrt(100*26-52) )
    
# Returns delta, the error in Stein's method.
def to_delta( covariances, length_genome, B_3, B_4 ):
    k = len(covariances)
    sigma_l = to_sigma_l( covariances, length_genome )
    D = 2*k-1
    radical = D*D*length_genome*B_3/sigma_l**3
    radical += 2*D**1.5*sqrt(7.0*length_genome*B_4)/sigma_l**2/sqrt(pi)
    delta = 2.0*(2.0/pi)**0.25*sqrt(radical)
    return delta
def _test_to_delta():
    B3 = B4 = 1.0
    radical = 7*7*100*B3/sqrt(2548)**3
    radical += 2*7**1.5*sqrt(7.0*100*B4)/sqrt(2548)**2/sqrt(pi)
    delta = 2.0*(2.0/pi)**0.25*sqrt(radical)
    assert isclose( to_delta( [2,3,4,5], 100, 1.0, 1.0 ), delta )

# Returns z_0.5_alpha, the normal quantile required in the central limit theorem.
#     Returns None unless delta < alpha0 to indicate the inequality is useless.
def z_half_alpha( alpha0, delta=0.0 ): # alpha0 is the desired size of the test
    if alpha0 <= delta:
        return None
    alpha = alpha0-delta
    z_half_alpha = norm.ppf(1.0-0.5*alpha)
    return z_half_alpha    
def _test_z_half_alpha():
    alpha0 = 0.15
    delta = 0.1
    assert isclose( z_half_alpha( alpha0, delta ), 1.959963984540054 )

def main(): 
    _test_to_sigma_square_x()
    _test_to_gamma_x()
    _test_to_sigma_l()
    _test_to_delta()
    _test_z_half_alpha()
    
if __name__ == "__main__":
    main()
