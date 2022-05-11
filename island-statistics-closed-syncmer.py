#!/usr/bin/env python
"""
Calculates confidence interval for sequence length from count of closed syncmers.
"""

# Patterned after:
#   https://github.com/medvedevgroup/mutation-rate-intervals

# Syncmer_Closed(k, s) # k-mer syncmers with s-codes

import sys
sys.path.append("./modules")

from argparse import ArgumentParser, RawTextHelpFormatter
from jls_syncmer_closed import Syncmer_Closed

def main(): 
    parser = getArguments()
    argument = parser.parse_args()
    check( argument )  

    k = argument.kmer_length
    s = argument.smer_length

    submer = Syncmer_Closed( k, s )
    p = submer.probability()

    h_output = {}
    h_output["islands:genome_fraction_of_starts"] = submer.first_passage_moment(0,k)*p
    h_output["islands:genome_fraction_of_bases"] = submer.first_passage_moment(1,k)*p
    h_output["islands:expected_size-biased_size"] = submer.first_passage_moment(2,k)*p
    output( h_output )    

# Outputs island statistics. 
def output( h_output ):
    print( *(h_output.keys()), sep='\t' )    
    print( *(h_output.values()), sep='\t' )

# Checks and fixes arguments if possible.    
def check( argument ):
    if argument.kmer_length <= 0:
        raise Exception( f'Error: kmer_length = {argument.kmer_length} <= 0.' )
    if argument.smer_length <= 0:
        raise Exception( f'Error: smer_length = {argument.smer_length} <= 0.' )
    if argument.kmer_length <= argument.smer_length:
        raise Exception( f'Error: smer_length <= kmer_length : {argument.smer_length} <= {argument.kmer_length}.' )
        
def getArguments():
    parser = ArgumentParser(description='Calculates island statistics for closed syncmers.\n',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument("-k", "--kmer_length", dest="kmer_length", type=int, required=True, 
                        help="KMER_LENGTH is the closed syncmer length k.", metavar="KMER_LENGTH")
    parser.add_argument("-s", "--smer_length", dest="smer_length", type=int, required=True,
                        help="SMER_LENGTH is the length s of s-mers within a closed syncmer.")
    return parser
    
if __name__ == "__main__":
    main()
