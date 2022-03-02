#!/usr/bin/env python
"""
Calculates confidence interval for mutation probability per letter from count of open syncmers, all and unmutated.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

import sys
sys.path.append("./modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from math import isclose
import jls_submer_clt_mgr
from jls_submer_to_interval_util import to_theta_interval
import jls_pattern_util

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    count_of_unmutated_submers = argument.number_of_unmutated_submers
    confidence = argument.confidence
    alpha_0 = 1.0 - confidence # for consistency with z-score
    length = argument.length
    
    k = argument.kmer_length
    patterns = argument.patterns
    nucleotide2probability = argument.nucleotide2probability

    try:
        h = jls_submer_clt_mgr.to_non_overlapping_pattern_Wilson_score_intervals_for_theta( count_of_submers, count_of_unmutated_submers, alpha_0, k, patterns, nucleotide2probability, length )
        interval = to_theta_interval( h )
    except:
        interval = [None, None]
    h_output = {}
    h_output["sig"] = confidence
    if length: 
        h_output["length"] = "actual"
    else:
        h_output["length"] = "estimated"
    h_output["ThetaLow"] = interval[0]
    h_output["ThetaHigh"] = interval[1]
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
    if argument.length is not None and argument.length <= 0:
        raise Exception( f'Error: length = {argument.length} <= 0.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates mutation probability per base from submer count of open syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-p", "--patterns", dest="patterns", type=str, required=True, 
                        help="PATTERNS is a set of patterns as IUPAC DNA strings, separated by commas (if 2 or more).", metavar="PATTERNS")
    parser.add_argument("-f", "--frequencies", dest="frequencies", type=str, default="0.25,0.25,0.25,0.25",
                        help="FREQUENCIES gives the probabilities for an independent letters model, in order, of acgt, defaulting to uniform frequencies.", metavar="FREQUENCIES")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of submers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-u", "--number_of_unmutated_submers", dest="number_of_unmutated_submers", type=int, required=True,
                        help="NUMBER_OF_UNMUTATED_SUBMERS is the count of unmutated submers within a sequence.", metavar="NUMBER_OF_UNMUTATED_SUBMERS")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE is the probability that the mutation probability per base lies within the interval [ThetaLow, ThetaHigh].", metavar="CONFIDENCE")
    parser.add_argument("-l", "--length", dest="length", type=float,  default=None, # ? use estimate for L ?
                        help="LENGTH (optional) inputs the known sequence length, if available.", metavar="LENGTH")
    return parser
    
if __name__ == "__main__":
    main()
