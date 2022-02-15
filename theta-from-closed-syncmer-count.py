#!/usr/bin/env python
"""
Calculates confidence interval for mutation probability per letter from count of closed syncmers, all and unmutated.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

from argparse import ArgumentParser, RawTextHelpFormatter
import jls_submer_clt_mgr
from jls_submer_to_interval_util import to_theta_interval

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  
    count_of_submers = argument.number_of_submers
    count_of_unmutated_syncmers = argument.number_of_unmutated_submers
    confidence = argument.confidence
    alpha_0 = 1.0 - confidence # for consistency with z-score
    length = argument.length
    k = argument.kmer_length
    s = argument.smer_length
    try:
        h = jls_submer_clt_mgr.to_syncmer_closed_Wilson_score_intervals_for_theta( count_of_submers, count_of_unmutated_syncmers, alpha_0, k, s, length )
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
        raise Exception( f'Error: smer_length <= kmer_length : {argument.smer_length} <= {argument.kmer_length}.' )
    if argument.confidence <= 0.0:
        raise Exception( f'Error: confidence = {argument.confidence} <= 0.0.' )
    if 1.0 <= argument.confidence:
        raise Exception( f'Error: 1.0 <= confidence = {argument.confidence}.' )
    if argument.length is not None and argument.length <= 0:
        raise Exception( f'Error: length = {argument.length} <= 0.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates mutation probability per base from submer count of closed syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the submer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within a syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-n", "--number_of_submers", dest="number_of_submers", type=int, required=True,
                        help="NUMBER_OF_SUBMERS is the count of submers within a sequence.", metavar="NUMBER_OF_SUBMERS")
    parser.add_argument("-u", "--number_of_unmutated_submers", dest="number_of_unmutated_submers", type=int, required=True,
                        help="NUMBER_OF_UNMUTATED_SUBMERS is the count of unmutated submers within a sequence.", metavar="NUMBER_OF_UNMUTATED_SUBMERS")
    parser.add_argument("-c", "--confidence", dest="confidence", type=float, default=0.95, # size of the confidence interval 
                        help="CONFIDENCE is the probability that the mutation probability per base lies within the interval [ThetaLow, ThetaHigh].", metavar="CONFIDENCE")
    parser.add_argument("-l", "--length", dest="length", type=float,  default=None, # ? use estimate for L ?
                        help="LENGTH (optional) uses an estimate of the sequence length, if available.", metavar="LENGTH")
    return parser
    
if __name__ == "__main__":
    main()
