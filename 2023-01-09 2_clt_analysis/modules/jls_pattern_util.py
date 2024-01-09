#!/usr/bin/env python
"""
Converts ambiguous nucleotides to regexes.
Internal representation is lower case, but interfaces can be specified as lower- or upper-case.
"""

from copy import deepcopy
import collections

# M.C. Frith, L. Noé, G. Kucherov (2020) 
# Minimally-overlapping words for sequence similarity search. 
# Bioinformatics 36: 5344-50.

# Returns regex corresponding to an ambiguous nucleotide pattern.
def to_regex(patterns): # Argument patterns may be a single string or a list of strings.
    regex = ''
    if isinstance(patterns, str):
        for letter in patterns: # single pattern
            unambiguous = iupac2nucleotides()[letter]
            if len(unambiguous) == 1:
                regex += unambiguous
            else:
                regex += '(' + '|'.join(iupac2nucleotides()[letter]) + ')'
        return '(' + regex + ')'
    elif isinstance(patterns, list): # list of patterns
        for pattern in patterns:
            if regex:
                regex += '|'
            regex += to_regex(pattern)
        return '(' + regex + ')'
    else:
        raise Exception(f'The patterns must be a IUPAC string or a list of IUPAC strings : {patterns}.')
def _test_to_regex():
    assert to_regex('cst') == '(c(c|g)t)'
    assert to_regex('ryy') == '((a|g)(c|t)(c|t))'
    assert to_regex(['a', 'cst']) == '((a)|(c(c|g)t))'
    assert to_regex(['a', 'cst']) == '((a)|(c(c|g)t))'
    # upper case
    assert to_regex('CST') == '(C(C|G)T)'
    assert to_regex('cSt') == '(c(C|G)t)'
    assert to_regex(['a', 'cSt']) == '((a)|(c(C|G)t))'
    
# Returns unambiguous DNA nucleotides as a string (length 4).
def nucleotides(): 
    return 'acgt'+'acgt'.upper()
# Returns IUPAC DNA nucleotides including ambiguous nucleotides as a dictionary (length 15).
def iupac2nucleotides(): # DNA nucleotides
    iupac2nucleotides = {
        'a':'a', 'c':'c', 'g':'g', 't':'t', 
        'r':'ag', 'y':'ct', 's':'cg', 'w':'at', 'k':'gt', 'm':'ac',
        'b':'cgt', 'd':'agt', 'h':'act', 'v':'acg', 
        'n':'acgt'
    }
    # Adds upper case letters to dictionary.
    lc = deepcopy(list(iupac2nucleotides.keys()))
    for n in lc: 
        iupac2nucleotides[n.upper()] = iupac2nucleotides[n].upper()
    return iupac2nucleotides
# ? Is string all DNA IUPAC ?
def is_dna_iupac(s): 
    IUPAC2NUCLEOTIDES = iupac2nucleotides()
    return all([letter in IUPAC2NUCLEOTIDES for letter in s])
def _test_is_dna_iupac():
    assert is_dna_iupac('a')
    assert is_dna_iupac('b')
    assert not is_dna_iupac('u')
    assert is_dna_iupac('abcdghtvrywskmn')
    assert not is_dna_iupac('abcdghtvrywskmnu')
    # upper case
    assert not is_dna_iupac('U')
    assert is_dna_iupac('ABCDGHTVRYWSKMN')
    assert is_dna_iupac('abcdGhtvrYwskmn')
# In an independent letters model, nucleotide2probability gives the probability of each of the 4 unambiguous nucleotide.
def to_probability(iupac, nucleotide2probability):
    nucleotides = iupac2nucleotides()[iupac]
    return sum([nucleotide2probability[nucleotide] for nucleotide in nucleotides])
def _test_to_probability():
    nucleotide2probability = {'a':0.1, 'c':0.2, 'g':0.3, 't':0.4}
    assert to_probability('a', nucleotide2probability) == 0.1
    assert to_probability('w', nucleotide2probability) == 0.5
    assert to_probability('n', nucleotide2probability) == 1.0
# Returns the string of lower-case unambiguous nucleotides common to IUPAC letters iupac0 & iupac.
def intersect(dna_iupac0, dna_iupac):
    NUCLEOTIDES = nucleotides()
    IUPAC2NUCLEOTIDES = iupac2nucleotides()
    intersect = ''
    for n in NUCLEOTIDES:
        if n.islower() and n in IUPAC2NUCLEOTIDES[dna_iupac.lower()] and n in IUPAC2NUCLEOTIDES[dna_iupac0.lower()]:
            intersect += n
    return intersect
def _test_intersect():
    assert intersect('a', 'b') == ''
    assert intersect('a', 'r') == 'a'
    assert intersect('a', 'y') == ''
    assert intersect('a', 'n') == 'a'
    assert intersect('d', 'b') == 'gt'
    assert intersect('r', 'y') == ''
    assert intersect('d', 'y') == 't'
    assert intersect('r', 'n') == 'ag'
    assert intersect('n', 'n') == 'acgt'
    # upper case
    assert intersect('N', 'n') == 'acgt'
    assert intersect('N', 'N') == 'acgt'
# Returns shortest postfix == prefix if pattern0[-j:] == pattern[i-j].
#   Returns '' otherwise.
#   Returns a non-empty string if pattern0 == pattern.
# The pattern0 & pattern are interpreted as IUPAC strings, 
#   so ANY common unambiguous nucleotide string in lower case matches both, e.g.,
#   postfix_eq_prefix('wS', 'ry', 2) = 'wS' because (lower case) 'ac' is common to both IUPAC strings.
def postfix_eq_prefix(pattern0, pattern, i=None): # i is the length of the prefix pattern[:i+1]
    length = min(len(pattern0), len(pattern))
    if i is None:
        for i in range(1,length+1):
            overlap = postfix_eq_prefix(pattern0, pattern, i)
            if overlap:
                return overlap
    else:
        assert i <= length
        if all([ intersect(pattern0[-j-1],pattern[i-j-1]) for j in range(i) ]):
            return pattern0[-i:]
    return ''
def _test_postfix_eq_prefix():
    # Returns a non-empty string if pattern0 == pattern.
    assert postfix_eq_prefix('ryy', 'ryy') == 'ryy'
    assert postfix_eq_prefix('abb', 'abb') == 'abb'
    assert postfix_eq_prefix('wba', 'wba') == 'a'
    assert postfix_eq_prefix('rYy', 'ryy') == 'rYy'
    assert postfix_eq_prefix('ryy', 'rYy') == 'ryy'
    assert postfix_eq_prefix('rYy', 'rYy') == 'rYy'
    assert postfix_eq_prefix('Wba', 'wba') == 'a'
    assert postfix_eq_prefix('wbA', 'wba') == 'A'
    assert postfix_eq_prefix('WBA', 'wba') == 'A'
    # _test_intersect(), with the return from postfix_eq_prefix
    assert postfix_eq_prefix('a', 'b') == ''
    assert postfix_eq_prefix('a', 'r') == 'a'
    assert postfix_eq_prefix('a', 'y') == ''
    assert postfix_eq_prefix('a', 'n') == 'a'
    assert postfix_eq_prefix('d', 'b') == 'd'
    assert postfix_eq_prefix('r', 'y') == ''
    assert postfix_eq_prefix('d', 'y') == 'd'
    assert postfix_eq_prefix('r', 'n') == 'r'
    assert postfix_eq_prefix('n', 'n') == 'n'
    # random testing
    assert postfix_eq_prefix('aa', 'bb') == ''
    assert postfix_eq_prefix('daa', 'yyy') == ''
    assert postfix_eq_prefix('aa', 'rr') == 'a'
    assert postfix_eq_prefix('ww', 'rr') == 'w'
    assert postfix_eq_prefix('aw', 'yy') == 'w' # 't' is common to 'w' and 'y'
    assert postfix_eq_prefix('aa', 'nn') == 'a'
    assert postfix_eq_prefix('aa', 'nnn') == 'a'
    assert postfix_eq_prefix('aaa', 'nn') == 'a'
    assert postfix_eq_prefix('rab', 'n') == 'b'
    assert postfix_eq_prefix('n', 'rab') == 'n'
    assert postfix_eq_prefix('da', 'bw') == 'da'
    assert postfix_eq_prefix('ada', 'bwn') == 'da'
    assert postfix_eq_prefix('rry', 'ryy') == 'ry'
    # similar non-overlapping words
    assert postfix_eq_prefix('ryy', 'ry') == ''
    assert postfix_eq_prefix('ry', 'ryy') == 'ry'
    # standard non-overlapping words
    assert postfix_eq_prefix('aab', 'abb') == 'ab'
    assert postfix_eq_prefix('aab', 'ab') == 'ab'
# Returns a non-empty overlap if some overlaps exist.
#   Returns '' otherwise.
#   Returns any replicated pattern in the patterns.
# Argument patterns may be a string or list of strings.
# The pattern0 & pattern are interpreted as IUPAC strings, 
#   as in postfix_eq_prefix.
def overlap(patterns):
    if isinstance(patterns, str):
        patterns = [patterns]
    for pattern, count in collections.Counter(patterns).items():
        if count > 1:
            return pattern
    for pattern0 in patterns:
        for pattern in patterns:
            postfix = postfix_eq_prefix(pattern0, pattern)
            if postfix and (pattern0 != pattern or postfix != pattern): 
                return postfix
    return ''
def _test_overlap():
    # Tests behavior on a single pattern.
    assert overlap('ryy') == ''
    assert overlap(['ryy']) == ''
    assert overlap('rYy') == ''
    assert overlap(['rYy']) == ''
    assert overlap('RYY') == ''
    assert overlap(['RYY']) == ''
    assert overlap('abb') == ''
    assert overlap(['abb']) == ''
    assert overlap('wba') == 'a'
    assert overlap(['wba']) == 'a'
    # Returns replicated patterns.
    assert overlap(['ryy', 'ryy']) == 'ryy'
    assert overlap(['rYy', 'ryy']) == 'rYy'
    assert overlap(['RYY', 'ryy']) == 'RYY'
    assert overlap(['abb', 'abb']) == 'abb'
    assert overlap(['wba', 'wba']) == 'wba'
    # _test_intersect(), with the return from postfix_eq_prefix
    assert overlap(['ryy', 'ry']) == 'ry'
    assert overlap(['ry', 'ryy']) == 'ry'
    assert overlap(['a', 'cg']) == '' # All are not overlapping and not self overlapping.
    assert overlap(['a', 'cst']) == '' # All are not overlapping and not self overlapping.
    assert overlap(['wss', 'sw']) == 's'
    assert overlap(['wsS', 'sw']) == 'S'
    assert overlap(['www', 'ssss']) == 'w' # self-overlap 'www'
    assert overlap(['Www', 'ssss']) == 'w'
    assert overlap(['wwW', 'ssss']) == 'W'
    assert overlap(['aaa', 'ttt', 'ssss']) == 'a'
 
def _test():
    _test_is_dna_iupac()
    _test_to_regex()
    _test_to_probability()
    _test_intersect()
    _test_postfix_eq_prefix()
    _test_overlap()
    
def main(): 
    _test()
    
if __name__ == "__main__":
    main()
