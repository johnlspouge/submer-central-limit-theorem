#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of open syncmers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes

from argparse import ArgumentParser, RawTextHelpFormatter
import jls_submer_clt_mgr
from jls_submer_to_interval_util import to_length_interval

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    confidence = argument.confidence
    alpha_0 = 1.0 - confidence # for consistency with z-score
    is_estimate = argument.is_estimate
    k = argument.kmer_length
    s = argument.smer_length
    t = argument.t_offset_of_smer
    try:
        h = jls_submer_clt_mgr.to_syncmer_open_Wilson_score_intervals_for_length( count_of_submers, alpha_0, k, s, t, is_estimate )
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
    if argument.kmer_length <= 0:
        raise Exception( f'Error: kmer_length = {argument.kmer_length} <= 0.' )
    if argument.smer_length <= 0:
        raise Exception( f'Error: smer_length = {argument.smer_length} <= 0.' )
    if argument.kmer_length <= argument.smer_length:
        raise Exception( f'Error: kmer_length <= smer_length : {argument.kmer_length} <= {argument.smer_length}.' )
    if argument.t_offset_of_smer < 0:
        raise Exception( f'Error: t_offset_of_smer = {argument.t_offset_of_smer} < 0.' )
    u = argument.kmer_length - argument.smer_length
    if u < argument.t_offset_of_smer:
        raise Exception( f'Error: kmer_length - smer_length < t_offset_of_smer : {u} < {argument.t_offset_of_smer}.' )
    if argument.confidence <= 0.0:
        raise Exception( f'Error: confidence = {argument.confidence} <= 0.0.' )
    if 1.0 <= argument.confidence:
        raise Exception( f'Error: 1.0 <= confidence = {argument.confidence}.' )
    if not argument.is_estimate:
        argument.is_estimate = None
        
def getArguments():
    parser = ArgumentParser(description='Calculates sequence length from submer count of open syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within an open syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-t", "--t_offset_of_smer", dest="t_offset_of_smer", type=int, required=True,
                        help="T_OFFSET_OF_SMER is the offset of the minimum s-mer within an open syncmer [0 <= t <= k-s].")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of open syncmers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE is the probability that the sequence length lies within the interval [LengthLow, LengthHigh].", metavar="CONFIDENCE")
    parser.add_argument("-e", "--is_estimate", dest="is_estimate", default=False, action="store_true", # ? use estimate for L in Stein's method ?
                        help="IS_ESTIMATE uses the estimate method with the flag and the qualitative method without it.")
    return parser
    
if __name__ == "__main__":
    main()
