#!/usr/bin/env python
"""
Mutation probability estimation by Stein's method and a central limit theorem for short-range dependency.
    For input from the abstract class Syncmer jls_syncmer.
"""
# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881
# Returns the syncmer covariances.

from enum import Enum
from math import isclose
from numpy import allclose
from scipy.stats import norm
from scipy.optimize import bisect
from jls_syncmer_parametrized import Syncmer_Parametrized # for tests
from jls_submer_clt_util import to_sigma_l # Uses to_sigma_square_x, to_gamma_x implicitly.

# Returns the point estimate, the mutation probability per letter.
def theta_mean( count_of_syncmers, count_of_unmutated_syncmers, k ):
    return 1.0 - (count_of_unmutated_syncmers/count_of_syncmers)**(1.0/k)

# Returns the standardized_variate_w W(theta, length_genome) 
#   in the central limit theorem about theta.
def standardized_variate_w(
    theta, # mutation probability per letter
    syncmer, # object from the Syncmer abstract class
    length_genome, # genome length for CLT, possibly a probabilistic bound
    count_of_syncmers, 
    count_of_unmutated_syncmers): 
    theta_bar = 1.0-theta
    k = syncmer.get_k()
    mutation_covariances = [ syncmer.mutation_covariance(i,theta) for i in range(k) ]
    try:
        sigma_l = to_sigma_l( mutation_covariances, length_genome )
    except:
        raise Exception('square root failed to produce sigma_l')
    standardized_variate_w = count_of_unmutated_syncmers-theta_bar**k*count_of_syncmers
    if sigma_l == 0.0:
        if standardized_variate_w < 0.0:
            return float('-inf')
        elif standardized_variate_w > 0.0:
            return float('inf')
        else:
            return 0.0
            #raise Exception('The standardized_variate_w subroutine yields 0/0.')
    standardized_variate_w /= sigma_l
    return standardized_variate_w

def _test_standardized_variate_w():
    theta = 0.0008388608
    (k,s,ts) = (6,2,[3])
    syncmer = Syncmer_Parametrized(k,s,ts)
    length_genome = 1000000
    count_of_syncmers = 200000
    count_of_unmutated_syncmers = 199000
    z = standardized_variate_w(theta, syncmer, length_genome, count_of_syncmers, count_of_unmutated_syncmers)
    z0 = 0.11571079067802262
    assert isclose(z,z0)

# Returns an incremented theta.
def _increment(theta):
    end_mesh = 1.0e-10
    mesh = 0.01
    factor = 2.0
    if theta == 0.0:
        return end_mesh
    if theta < mesh:
        return theta*factor
    elif theta < 1.0-mesh:
        return theta+mesh
    elif theta < 1.0-end_mesh:
        return ((factor-1.0)*theta+1.0)/factor     
    else:
        theta = 1.0
    return theta

def _test_increment():
    theta = 0.0
    while theta != 1.0:
        theta0 = theta 
        theta = _increment(theta0)
        assert 0.0 <= theta <= 1.0
        assert 0.0 <= theta-theta0 <= 0.01 + 1.0e-06

# Returns the Wilson score intervals for approximation for theta, the mutation probability per letter.
def _theta_zero( theta_minus, theta_plus ):
    THETA_0 = 0.0
    THETA_1 = 1.0
    theta_zero = theta_minus+theta_plus
    theta_zero.sort()
    if theta_zero[0] == THETA_0:
        theta_zero.pop(0)
    else:
        theta_zero.insert(0, THETA_0)
    if theta_zero[-1] == THETA_1:
        theta_zero.pop()
    else:
        theta_zero.append(THETA_1)
    return theta_zero

# Returns the three Wilson score intervals for approximation for theta, the mutation probability per letter.
def theta_all(count_of_syncmers, count_of_unmutated_syncmers, z_half_epsilon, syncmer, length_genome ):
    theta_minus = []
    theta_zero  = []
    theta_plus  = []
    class Interval(Enum):
        MINUS = -1
        ZERO  =  0
        PLUS  =  1
    interval2list = {
        Interval.MINUS:theta_minus, 
        Interval.ZERO:theta_zero, 
        Interval.PLUS:theta_plus
    }
    def _interval(w):
        if w < -z_half_epsilon:
            return Interval.MINUS
        elif w < z_half_epsilon:
            return Interval.ZERO
        else:
            return Interval.PLUS
    def _update(cross0, interval0, theta0, interval, theta):
        if interval0 == Interval.MINUS or interval == Interval.MINUS: # Prioritizes interval0.
            value = -z_half_epsilon
        elif interval0 == Interval.PLUS or interval == Interval.PLUS: # Prioritizes interval0.
            value = z_half_epsilon
        else:
            raise Exception('The mesh is not fine enough to localize changes.')
        def _zero(theta):
            w = standardized_variate_w(theta, syncmer, length_genome, count_of_syncmers, count_of_unmutated_syncmers)
            return w-value
        cross = bisect( _zero, theta0, theta )
        if interval0 != Interval.ZERO: # Completes Interval.MINUS or Interval.PLUS.
            interval2list[interval0].extend([cross0, cross])
            if interval != Interval.ZERO:
                cross0 = cross
                if interval == Interval.MINUS:
                    value = -z_half_epsilon
                elif interval == Interval.PLUS:
                    value = z_half_epsilon
                else:
                    raise Exception('The mesh is not fine enough to localize changes.')
                cross = bisect( _zero, cross0, theta )
        cross0 = cross
        interval0 = interval
        return interval0, cross0
    interval0 = None # current interval
    cross0 = 0.0 # previous endpoint
    theta0 = 0.0 # current theta
    while theta0 != 1.0:
        theta = _increment( theta0 )
        w = standardized_variate_w(theta, syncmer, length_genome, count_of_syncmers, count_of_unmutated_syncmers)
        interval = _interval(w)
        if interval0 is None: # Initializes Interval, because theta=0.0 may be singular.
            interval0 = interval
        if interval0 != interval: # The function has crossed at least one boundary.
            interval0, cross0 = _update(cross0, interval0, theta0, interval, theta)
        theta0 = theta
    interval2list[interval0].extend([cross0, 1.0])
    return _theta_zero(theta_minus, theta_plus)

def _test_theta_all():
    theta_mean0s = [
        0.0008350747680885284,
        0.008512444610847103, # alpha = 0.05
        0.10910128185966073,
        0.31870793094203864
    ]
    theta_all0s = [
        [
            [0.0007734865796875001, 0.0009015489609375002],
            [0.0003947706015625, 0.0017629700265624998],
            [1.3862220312500001e-05, 0.0423956341611553],
        ],[  # alpha = 0.05
            [0.0083097677921875, 0.008719889891191962],
            [0.0066882563265625, 0.010812503117187499],
            [0.0010934997265624996, 0.05574933284552357],
        ],[
            [0.1083230310248542, 0.10988257118642374],
            [0.10146205195350158, 0.11704388602761136],
            [0.0500725099480111, 0.1944894128421653],
        ],[
            [0.3170349083992985, 0.3203777585756129],
            [0.30185928522966404, 0.3352391716966739],
            [0.16281319724520252, 0.4566894464114083],
        ]
    ]
    alpha = 0.05
    z_half_epsilon = norm.ppf(1.0-0.5*alpha)
    (k,s,ts) = (6,2,[3])
    syncmer = Syncmer_Parametrized(k,s,ts)
    fractions = [0.005, 0.05, 0.5, 0.9]
    length_genomes = [1000000, 10000, 100]
    for i in range(len(fractions)):
        for j in range(len(length_genomes)):
            fraction = fractions[i]
            length_genome = length_genomes[j]
            count_of_syncmers = length_genome/(k-s+1) # 200000
            count_of_unmutated_syncmers = count_of_syncmers-fraction*count_of_syncmers # 199000
            theta_all_ = theta_all( count_of_syncmers, count_of_unmutated_syncmers, z_half_epsilon, syncmer, length_genome )
            theta_mean_ = theta_mean( count_of_syncmers, count_of_unmutated_syncmers, k )
        isclose(theta_mean_, theta_mean0s[j])
        allclose(theta_all_,theta_all0s[i][j])

def main(): 
    _test_standardized_variate_w()
    _test_increment()
    _test_theta_all()
    
if __name__ == "__main__":
    main()
