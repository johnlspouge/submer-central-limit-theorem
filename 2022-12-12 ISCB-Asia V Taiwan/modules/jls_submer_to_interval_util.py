#!/usr/bin/env python
"""
Routines for converting hashes containing confidence intervals into a single interval (if possible).
"""

# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881

from math import isinf

# Returns confidence interval if hash of score intervals is typical.
#   Returns [None, None] otherwise.
def to_length_interval( confidence_length_hash ):
    if _is_typical_confidence_length_hash( confidence_length_hash ):
        return confidence_length_hash["confidence"]
    return [None, None]
        
# Returns confidence interval if hash of score intervals is typical.
#   Returns [None, None] otherwise.
def to_theta_interval( confidence_theta_hash ):
    if _is_typical_confidence_theta_hash( confidence_theta_hash ):
        return confidence_theta_hash["confidence"]
    return [None, None]
        
# Returns keys in typical hash of score intervals.
def _confidence_hash_keys():
    return ("left", "confidence", "right")

# Returns True if h is a hash containing keys from _confidence_hash_keys(), each a list of length 2.
def _is_typical_confidence_hash( h ):
    KEYS = _confidence_hash_keys()
    for key in KEYS:
        if key not in h or len( h[key] ) != 2 or h[key][1] <= h[key][0]:
            return False
    return True
    
# Returns True if h is a typical hash for score intervals of length.
def _is_typical_confidence_length_hash( h ):
    if not _is_typical_confidence_hash( h ):
        return False
    x = h["left"][1]
    if not isinf(x) or x < 0:
        return False
    if h["left"][0] != h["confidence"][1]:
        return False
    if h["confidence"][0] != h["right"][1]:
        return False
    return True
    
def _test_is_typical_confidence_length_hash():
    h = {'left': [2, float('inf')], 'confidence': [1, 2], 'right': [0, 1]}
    assert _is_typical_confidence_length_hash( h )
    # missing key
    h = {'confidence': [1, 2], 'right': [0, 1]}
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [float('inf'), 2], 'right': [0, 1]}
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [float('inf'), 2], 'confidence': [1, 2]}
    assert not _is_typical_confidence_length_hash( h )
    # interval not ordered
    h = {'left': [float('inf'), 2], 'confidence': [1, 2], 'right': [0, 1]}
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [2, float('inf')], 'confidence': [2, 1], 'right': [0, 1]}
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [2, float('inf')], 'confidence': [1, 2], 'right': [1, 0]}
    assert not _is_typical_confidence_length_hash( h )
    # intervals in incorrect order
    h = {'left': [2, 3], 'confidence': [1, 2], 'right': [0, 1]}
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [0, 1], 'confidence': [1, 2], 'right': [2, float('inf')] }
    assert not _is_typical_confidence_length_hash( h )
    h = {'left': [2, 3], 'confidence': [1, 2], 'right': [2, float('inf')] }
    assert not _is_typical_confidence_length_hash( h )
    
# Returns True if h is a typical hash for score intervals of theta.
def _is_typical_confidence_theta_hash( h ):
    if not _is_typical_confidence_hash( h ):
        return False
    if h["left"][0] != 0.0:
        return False
    if h["left"][1] != h["confidence"][0]:
        return False
    if h["confidence"][1] != h["right"][0]:
        return False
    if h["right"][1] != 1.0:
        return False
    return True
    
def _test_is_typical_confidence_theta_hash():
    h = {'left': [0.0, 0.1], 'confidence': [0.1, 0.9], 'right': [0.9, 1.0]}
    assert _is_typical_confidence_theta_hash( h )
    # missing key
    h = {'confidence': [0.1, 0.9], 'right': [0.9, 1.0]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.0, 0.1], 'right': [0.9, 1.0]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.0, 0.1], 'confidence': [0.1, 0.9]}
    assert not _is_typical_confidence_theta_hash( h )
    # interval not ordered
    h = {'left': [0.1, 0.0], 'confidence': [0.1, 0.9], 'right': [0.9, 1.0]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.0, 0.1], 'confidence': [0.9, 0.1], 'right': [0.9, 1.0]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.0, 0.1], 'confidence': [0.1, 0.9], 'right': [1.0, 0.9]}
    assert not _is_typical_confidence_theta_hash( h )
    # intervals in incorrect order
    h = {'left': [0.1, 0.9], 'confidence': [0.0, 0.1], 'right': [0.9, 1.0]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.0, 0.1], 'confidence': [0.9, 1.0], 'right': [0.1, 0.9]}
    assert not _is_typical_confidence_theta_hash( h )
    h = {'left': [0.9, 1.0], 'confidence': [0.1, 0.9], 'right': [0.0, 0.1]}
    assert not _is_typical_confidence_theta_hash( h )
    
def main(): 
    _test_is_typical_confidence_length_hash()
    _test_is_typical_confidence_theta_hash()
    
if __name__ == "__main__":
    main()
