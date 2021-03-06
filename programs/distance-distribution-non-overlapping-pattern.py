#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of open syncmers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes and t-offset

import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from math import isclose
from jls_submer import Submer
import jls_pattern_util
from jls_non_overlapping_pattern import Non_Overlapping_Pattern

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  

    k = argument.kmer_length
    patterns = argument.patterns
    nucleotide2probability = argument.nucleotide2probability
    
    m = argument.max_distance
    is_test_probabilities = argument.is_test_probabilities
    
    submer = Non_Overlapping_Pattern( k, patterns, nucleotide2probability )
    p = []
    for i in range(0,m+1):    
        p.append(submer.first_passage_probability(i))
    if is_test_probabilities:
        p = Submer.to_test_probabilities(submer.probability(), p)
    output( p )    

# Outputs array[0...m]. 
def output( p ):
    print( *p, sep='\t' )

# Check and fixes arguments if possible.    
def check( argument ):
    LEN_NUCLEOTIDES = 4
    if argument.kmer_length <= 0:
        raise Exception( f'Error: kmer_length = {argument.kmer_length} <= 0.' )
    # argument.patterns
    patterns = argument.patterns.split(',')
    for pattern in patterns:
        if not jls_pattern_util.is_dna_iupac(pattern):
            raise Exception( f'Error: {pattern} is not a DNA IUPAC string.' )
    argument.patterns = patterns
    # argument.frequencies
    frequencies = argument.frequencies.split(',')
    if len(frequencies) != LEN_NUCLEOTIDES:
        raise Exception( f'Error: Each of a,c,g,t (in that order) should be given a frequency : {argument.frequencies}.' )
    try:
        frequencies = list(map(float, frequencies))
    except:
        raise Exception( f'Error: Each of the {LEN_NUCLEOTIDES} frequencies should be a float: {argument.frequencies}.' )
    if not isclose(sum(frequencies), 1.0):
        raise Exception( f'Error: The {LEN_NUCLEOTIDES} frequencies should sum to 1.0: {argument.frequencies} and sum = {sum(frequencies)}.' )
    nucleotides = 'acgt'  
    nucleotide2probability = {}
    for i in range(LEN_NUCLEOTIDES):
        nucleotide2probability[nucleotides[i]] = frequencies[i]
    argument.nucleotide2probability = nucleotide2probability
    # remaining arguments
    if argument.max_distance <= 0:
        raise Exception( f'Error: max_distance <= 0 : {argument.max_distance} <= 0.' )
    if not argument.is_test_probabilities:
        argument.is_test_probabilities = None

def getArguments():
    parser = ArgumentParser(description='Calculates distance distribution for non-overlapping patterns.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the submer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-p", "--patterns", dest="patterns", type=str, required=True, 
                        help="PATTERNS is a set of patterns as IUPAC DNA strings, separated by commas (if 2 or more).", metavar="PATTERNS")
    parser.add_argument("-f", "--frequencies", dest="frequencies", type=str, default="0.25,0.25,0.25,0.25",
                        help="FREQUENCIES gives the probabilities for an independent letters model, in order, of acgt, defaulting to uniform frequencies.", metavar="FREQUENCIES")
    parser.add_argument("-m", "--max_distance", dest="max_distance", type=int, required=True,
                        help="MAX_DISTANCE is the maximum distance of interest between minimally overlapping k-mers.")
    parser.add_argument("-y", "--is_test_probabilities", dest="is_test_probabilities", default=False, action="store_true", # ? alpha-test probabilities ?
                        help="IS_TEST_PROBABILITIES calculates Shaw & Yu alpha-test probabilities with the flag and first-passage probabilities without it.")
    return parser
    
if __name__ == "__main__":
    main()
