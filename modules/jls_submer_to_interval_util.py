#!/usr/bin/env python
"""
Routines for converting hashes containing confidence intervals into a single interval (if possible).
"""

# A. Blanca et al. (2021) 
# The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches. 
# DOI: 10.1101/2021.01.15.426881

# Returns confidence interval if hash of score intervals is typical.
#   Returns [None, None] otherwise.
def to_length_interval( confidence_length_hash ):
    if _is_typical_confidence_length_hash(confidence_length_hash):
        return confidence_length_hash
    return [None, None]
        
# Returns confidence interval if hash of score intervals is typical.
#   Returns [None, None] otherwise.
def to_theta_interval( confidence_theta_hash ):
    if _is_typical_confidence_theta_hash( confidence_theta_hash ):
        return confidence_theta_hash
    return [None, None]
        
# Returns True if h is a typical hash for score intervals of length.
def _is_typical_confidence_length_hash( h ):
    if len(h) != 2:
        return False
    if not isinstance(h[0],(float,int)) or not isinstance(h[1],(float,int)):
        return False
    if not 0.0 <= h[0] <= h[1]:
        return False
    return True
    
def _test_is_typical_confidence_length_hash():
    # interval not ordered
    h = [2, 1]
    assert not _is_typical_confidence_length_hash( h )
    
# Returns True if h is a typical hash for score intervals of theta.
def _is_typical_confidence_theta_hash( h ):
    if len(h) != 2:
        return False
    if not isinstance(h[0],(float,int)) or not isinstance(h[1],(float,int)):
        return False
    if not 0.0 <= h[0] <= h[1] <= 1.0:
        return False
    return True
    
def _test_is_typical_confidence_theta_hash():
    # interval not ordered
    h = [0.9, 0.1]
    assert not _is_typical_confidence_theta_hash( h )
    # typical interval
    h = [0.0, 0.1]
    assert _is_typical_confidence_theta_hash( h )
    h = [0.1, 0.9]
    assert _is_typical_confidence_theta_hash( h )
    h = [0.9, 1.0]
    assert _is_typical_confidence_theta_hash( h )
    
def main(): 
    _test_is_typical_confidence_length_hash()
    _test_is_typical_confidence_theta_hash()
    
if __name__ == "__main__":
    main()
