#!/usr/bin/env python
"""
Calculates hypothesis test p-values for mutation probability per letter from count of parametrized syncmers, all and unmutated.
"""
# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals
# Syncmer_Synchronized(k, s, ts) # k-mer syncmers with s-codes and parameters ts
import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
import scipy.stats as st
from math import isclose

from jls_submer_clt_theta import standardized_variate_w
from jls_non_overlapping_pattern import Non_Overlapping_Pattern
import jls_pattern_util

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    count_of_unmutated_submers = argument.number_of_unmutated_submers
    theta0 = argument.theta0
    length = argument.length
    
    k = argument.kmer_length
    patterns = argument.patterns
    nucleotide2probability = argument.nucleotide2probability
    non_overlapping_pattern = Non_Overlapping_Pattern(k, patterns, nucleotide2probability)
    
    h_output = {}
    if length: 
        h_output["length"] = "actual"
    else:
        h_output["length"] = "estimated"
        length = count_of_submers/non_overlapping_pattern.probability()
    try:
        w = standardized_variate_w(theta0, non_overlapping_pattern, length, count_of_submers, count_of_unmutated_submers)
        q = st.norm.cdf(w)
        h_output["p_1-sided_left"] = q
        h_output["p_1-sided_right"] = 1.0-q
        q0 = min(q,1.0-q)
        h_output["p_2-sided"] = 2.0*q0
    except:
        h_output["p_1-sided_left"] = h_output["p_1-sided_right"] = h_output["p_2-sided"] = None
    output( h_output )    

# Outputs confidence interval of genomic length. 
#   Outputs None as the p-values if the computation has an abnormal result.
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
    if argument.number_of_submers <= 0:
        raise ValueError( f'Error: length = {argument.number_of_submers} <= 0.' )
    if argument.number_of_unmutated_submers is not None and argument.number_of_unmutated_submers <= 0:
        raise ValueError( f'Error: length = {argument.length} <= 0.' )
    if not 0.0 <= argument.theta0 <= 1.0:
        raise ValueError( f'Error: 0.0 <= theta0 = {argument.theta0} <= 1.0.' )
    if argument.length is not None and argument.length <= 0:
        raise ValueError( f'Error: length = {argument.length} <= 0.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates hypothesis test p-values for mutation probability per base from submer count of non-overlapping patterns.\n',
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
    parser.add_argument("-0", "--theta0", dest="theta0", type=float, # theta0, the value of theta for hypothesis testing 
                        help="THETA0 is the mutation probability per base for hypothesis testing.", metavar="THETA0")
    parser.add_argument("-l", "--length", dest="length", type=float,  default=None, # ? use estimate for L ?
                        help="LENGTH (optional) inputs the known sequence length, if available.", metavar="LENGTH")
    return parser
    
if __name__ == "__main__":
    main()
