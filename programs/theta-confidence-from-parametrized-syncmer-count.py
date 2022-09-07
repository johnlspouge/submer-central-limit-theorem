#!/usr/bin/env python
"""
Calculates confidence interval for mutation probability per letter from count of parametrized syncmers, all and unmutated.
"""
# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals
# Syncmer_Synchronized(k, s, ts) # k-mer syncmers with s-codes and parameters ts
import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from jls_submer_clt_mgr import to_syncmer_parametrized_Wilson_score_intervals_for_theta
from jls_submer_to_interval_util import to_theta_interval

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_syncmers = argument.number_of_submers
    count_of_unmutated_syncmers = argument.number_of_unmutated_submers
    confidence = argument.confidence
    alpha_0 = 1.0 - confidence # for consistency with z-score
    length = argument.length
    k = argument.kmer_length
    s = argument.smer_length
    ts = argument.t_offsets_of_smer
    eps = argument.eps
    try:
        h = to_syncmer_parametrized_Wilson_score_intervals_for_theta( count_of_syncmers, count_of_unmutated_syncmers, alpha_0, k, s, ts, eps, length )
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
                        help="KMER_LENGTH is the submer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within a syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-t", "--t_offsets_of_smer", dest="t_offsets_of_smer", type=int, nargs='+', required=True,
                        help="T_OFFSETS_OF_SMER is the offset(s) of the minimum s-mer within a paramtrized syncmer.", metavar="T_OFFSETS_OF_SMER")
    parser.add_argument("-e", "--eps", dest="eps", type=float, default=0.0,
                        help="EPS is the rejection probability when downsampling paramtrized syncmers (EPS=0.0 for no downsampling).", metavar="EPS")
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
