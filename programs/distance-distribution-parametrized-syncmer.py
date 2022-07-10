#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of open syncmers.
"""
# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals
# Syncmer_Synchronized(k, s, ts) # k-mer syncmers with s-codes and parameters ts
import sys
sys.path.insert(0,"../modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from jls_submer import Submer
from jls_syncmer_parametrized import Syncmer_Parametrized

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  

    k = argument.kmer_length
    s = argument.smer_length
    ts = argument.t_offsets_of_smer
    m = argument.max_distance
    is_test_probabilities = argument.is_test_probabilities
    
    submer = Syncmer_Parametrized( k, s, ts )
    p = []
    for i in range(0,m+1):    
        p.append(submer.first_passage_probability(i))
    if is_test_probabilities:
        p = Submer.to_test_probabilities(submer.probability(), p)
    output( p )    

# Outputs array[0...m]. 
def output( p ):
    print( *p, sep='\t' )

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
    if argument.max_distance <= 0:
        raise Exception( f'Error: max_distance <= 0 : {argument.max_distance} <= 0.' )
    if not argument.is_test_probabilities:
        argument.is_test_probabilities = None
        
def getArguments():
    parser = ArgumentParser(description='Calculates distance distribution for open syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the open syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within an open syncmer.", metavar="SMER_LENGTH")
    parser.add_argument("-t", "--t_offsets_of_smer", dest="t_offsets_of_smer", type=int, nargs='+', required=True,
                        help="T_OFFSETS_OF_SMER is the offset(s) of the minimum s-mer within a paramtrized syncmer.", metavar="T_OFFSETS_OF_SMER")
    parser.add_argument("-m", "--max_distance", dest="max_distance", type=int, required=True,
                        help="MAX_DISTANCE is the maximum distance of interest between open syncmers.")
    parser.add_argument("-y", "--is_test_probabilities", dest="is_test_probabilities", default=False, action="store_true", # ? alpha-test probabilities ?
                        help="IS_TEST_PROBABILITIES calculates Shaw & Yu alpha-test probabilities with the flag and first-passage probabilities without it.")
    return parser
    
if __name__ == "__main__":
    main()
