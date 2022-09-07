#!/usr/bin/env python
"""
Calculates hypothesis test p-values for sequence length from count of parametrized syncmers.
"""
# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals
# Syncmer_Synchronized(k, s, ts) # k-mer syncmers with s-codes and parameters ts
import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
import scipy.stats as st

from jls_submer_clt_genome_length import standardized_variate_w
from jls_minimizer import Minimizer

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    length0 = argument.length0
    
    k = argument.kmer_length
    w = argument.window_length
    minimizer = Minimizer(w,k)
    
    h_output = {}
    try:
        w = standardized_variate_w(length0, minimizer, count_of_submers)
        q = st.norm.cdf(w)
        h_output["p_1-sided_left"] = q
        h_output["p_1-sided_right"] = 1.0-q
        q0 = min(q,1.0-q)
        h_output["p_2-sided"] = 2.0*q0
    except:
        h_output["p_1-sided_left"] = h_output["p_1-sided_right"] = h_output["p_2-sided"] = None
    output( h_output )    

# Outputs confidence interval of genomic length. 
#   Outputs [None, None] as the interval if the computation has an abnormal result.
def output( h_output ):
    print( *(h_output.keys()), sep='\t' )    
    print( *(h_output.values()), sep='\t' )

# Check and fixes arguments if possible.    
def check( argument ):
    if argument.kmer_length <= 0:
        raise Exception( f'Error: kmer_length = {argument.kmer_length} < 0.' )
    if argument.window_length <= 1:
        raise Exception( f'Error: window_length = {argument.window_length} <= 1.' )
    # remaining arguments
    if argument.number_of_submers <= 0:
        raise ValueError( f'Error: length = {argument.number_of_submers} <= 0.' )
    if argument.length0 <= 0:
        raise ValueError( f'Error: length = {argument.length} <= 0.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates hypothesis test p-values for sequence length from submer count of minimizers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the minimizer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-w", "--window_length", dest="window_length", type=int, required=True,
                        help="WINDOW_LENGTH is the number of k-mers in each minimizer window.", metavar="WINDOW_LENGTH")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of submers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-0", "--length0", dest="length0", type=float, # # length0, the value of length for hypothesis testing 
                        help="LENGTH0 is the sequence length for hypothesis testing.", metavar="LENGTH0")
    return parser
    
if __name__ == "__main__":
    main()
