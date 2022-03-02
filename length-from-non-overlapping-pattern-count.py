#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of open syncmers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes

import sys
sys.path.append("./modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from math import isclose
import jls_submer_clt_mgr
from jls_submer_to_interval_util import to_length_interval
import jls_pattern_util

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument ) 
    
    count_of_submers = argument.number_of_submers
    confidence = argument.confidence
    alpha_0 = 1.0 - confidence # for consistency with z-score
    is_estimate = argument.is_estimate

    k = argument.kmer_length
    patterns = argument.patterns
    nucleotide2probability = argument.nucleotide2probability

    try:
        h = jls_submer_clt_mgr.to_non_overlapping_pattern_Wilson_score_intervals_for_length( count_of_submers, alpha_0, k, patterns, nucleotide2probability, is_estimate )
        interval = to_length_interval( h )
    except:
        interval = [None, None]
    h_output = {}
    h_output["sig"] = confidence
    h_output["method"] = "qualitative"
    if is_estimate: 
        h_output["method"] = "estimate"
    h_output["LengthLow"] = interval[0]
    h_output["LengthHigh"] = interval[1]
    output( h_output )    

# Outputs confidence interval of genomic length. 
#   Outputs [None, None] as the interval if the computation has an abnormal result.
def output( h_output ):
    print( *(h_output.keys()), sep='\t' )    
    print( *(h_output.values()), sep='\t' )

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
    if argument.confidence <= 0.0:
        raise Exception( f'Error: confidence = {argument.confidence} <= 0.0.' )
    if 1.0 <= argument.confidence:
        raise Exception( f'Error: 1.0 <= confidence = {argument.confidence}.' )
    if not argument.is_estimate:
        argument.is_estimate = None
        
def getArguments():
    parser = ArgumentParser(description='Calculates sequence length from submer count of non-overlapping patterns.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-p", "--patterns", dest="patterns", type=str, required=True, 
                        help="PATTERNS is a set of patterns as IUPAC DNA strings, separated by commas (if 2 or more).", metavar="PATTERNS")
    parser.add_argument("-f", "--frequencies", dest="frequencies", type=str, default="0.25,0.25,0.25,0.25",
                        help="FREQUENCIES gives the probabilities for an independent letters model, in order, of acgt, defaulting to uniform frequencies.", metavar="FREQUENCIES")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of open syncmers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE is the probability that the sequence length lies within the interval [LengthLow, LengthHigh].", metavar="CONFIDENCE")
    parser.add_argument("-e", "--is_estimate", dest="is_estimate", default=False, action="store_true", # ? use estimate for L in Stein's method ?
                        help="IS_ESTIMATE uses the estimate method with the flag and the qualitative method without it.")
    return parser
    
if __name__ == "__main__":
    main()
