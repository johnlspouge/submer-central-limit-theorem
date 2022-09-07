#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of minimizers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Minimizer(w, k) # (w,k)-minimizers = w-windows of k-mers

import sys
sys.path.insert(0,"../modules")

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
    w = argument.window_length
    try:
        h = jls_submer_clt_mgr.to_minimizer_Wilson_score_intervals_for_length( count_of_submers, alpha_0, w, k, is_estimate )
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
        raise Exception( f'Error: kmer_length = {argument.kmer_length} < 0.' )
    if argument.window_length <= 1:
        raise Exception( f'Error: window_length = {argument.window_length} <= 1.' )
    # remaining arguments
    if argument.number_of_submers <= 0:
        raise ValueError( f'Error: length = {argument.number_of_submers} <= 0.' )
    if argument.confidence <= 0.0:
        raise Exception( f'Error: confidence = {argument.confidence} <= 0.0.' )
    if 1.0 <= argument.confidence:
        raise Exception( f'Error: 1.0 <= confidence = {argument.confidence}.' )
    if not argument.is_estimate:
        argument.is_estimate = None
        
def getArguments():
    parser = ArgumentParser(description='Calculates sequence length from the count of minimizers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the minimizer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-w", "--window_length", dest="window_length", type=int, required=True,
                        help="WINDOW_LENGTH is the number of k-mers in each minimizer window.", metavar="WINDOW_LENGTH")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of minimizers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE is the probability that the sequence length lies within the interval [LengthLow, LengthHigh].", metavar="CONFIDENCE")
    parser.add_argument("-e", "--is_estimate", dest="is_estimate", default=False, action="store_true", # ? use estimate for L in Stein's method ?
                        help="IS_ESTIMATE uses the estimate method with the flag and the qualitative method without it.")
    return parser
    
if __name__ == "__main__":
    main()
