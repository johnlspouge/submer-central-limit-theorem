#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of open syncmers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Syncmer_Open(k, s, t) # k-mer syncmers with s-codes and t-offset

import sys
sys.path.append("./modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from jls_syncmer_open import Syncmer_Open

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  

    k = argument.kmer_length
    s = argument.smer_length
    t = argument.t_offset_of_smer
    
    submer = Syncmer_Open( k, s, t )
    p = submer.probability()

    h_output = {}
    h_output["islands:genome_fraction_of_starts"] = submer.first_passage_moment(0,k)*p
    h_output["islands:genome_fraction_of_bases"] = submer.first_passage_moment(1,k)*p
    h_output["islands:expected_size-biased_size"] = submer.first_passage_moment(2,k)*p
    output( h_output )    

# Outputs confidence interval of genomic length. 
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
        
def getArguments():
    parser = ArgumentParser(description='Calculates distance distribution for open syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within an open syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-t", "--t_offset_of_smer", dest="t_offset_of_smer", type=int, required=True,
                        help="T_OFFSET_OF_SMER is the offset of the minimum s-mer within an open syncmer [0 <= t <= k-s].")
    return parser
    
if __name__ == "__main__":
    main()
