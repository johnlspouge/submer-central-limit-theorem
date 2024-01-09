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

from jls_syncmer_parametrized import Syncmer_Parametrized
from jls_submer_clt_genome_length import standardized_variate_w

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    length0 = argument.length0
    k = argument.kmer_length
    s = argument.smer_length
    ts = argument.t_offsets_of_smer
    eps = argument.eps
    syncmer = Syncmer_Parametrized(k, s, ts, eps)

    h_output = {}
    try:
        w = standardized_variate_w(length0, syncmer, count_of_submers)
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
        raise Exception( f'Error: kmer_length = {argument.kmer_length} <= 0.' )
    if argument.smer_length <= 0:
        raise Exception( f'Error: smer_length = {argument.smer_length} <= 0.' )
    if argument.kmer_length <= argument.smer_length:
        raise Exception( f'Error: kmer_length <= smer_length : {argument.kmer_length} <= {argument.smer_length}.' )
    u = argument.kmer_length - argument.smer_length
    for t in argument.t_offsets_of_smer:
        if t < 0:
            raise Exception( f'Error: t_offsets_of_smer = {argument.t_offsets_of_smer} < 0.' )
        if u < t:
            raise Exception( f'Error: kmer_length - smer_length < t_offset_of_smer : {u} < {argument.t_offsets_of_smer}.' )
    if not 0.0 <= argument.eps < 1.0:
        raise Exception( f'Error: the downsampling probability must satisfy 0.0 <= eps = {argument.eps} < 1.0.' )
    # remaining arguments
    if argument.number_of_submers <= 0:
        raise ValueError( f'Error: length = {argument.number_of_submers} <= 0.' )
    if argument.length0 <= 0:
        raise ValueError( f'Error: length = {argument.length} <= 0.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates hypothesis test p-values for sequence length from submer count of parametrized syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within an open syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-t", "--t_offsets_of_smer", dest="t_offsets_of_smer", type=int, nargs='+', required=True,
                        help="T_OFFSETS_OF_SMER is the offset(s) of the minimum s-mer within a paramtrized syncmer.", metavar="T_OFFSETS_OF_SMER")
    parser.add_argument("-e", "--eps", dest="eps", type=float, default=0.0,
                        help="EPS is the rejection probability when downsampling paramtrized syncmers (EPS=0.0 for no downsampling).", metavar="EPS")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of open syncmers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-0", "--length0", dest="length0", type=float, required=True, # # length0, the value of length for hypothesis testing 
                        help="LENGTH0 is the sequence length for hypothesis testing.", metavar="LENGTH0")
    return parser
    
if __name__ == "__main__":
    main()
