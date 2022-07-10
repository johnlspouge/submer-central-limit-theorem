# submer-central-limit-theorem

## The central limit theorems  (**CLT**s)

The code was inspired by the preprint of

A. Blanca et al. (2021)<br/> 
The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches.<br/> 
Journal of Computational Biology 29:155-168
DOI: 10.1089/cmb.2021.0431

which used Stein's method to develop CLTs for the k-mer sketch of sequences. Under reasonable assumptions, the k-mer sketch can determine sequence length exactly because the it lists all k-mers in the sequence. In contrast, submers sample the k-mers in a sequence. Thus, the submer count does not determine the sequence length exactly, but yields probabilistic bounds on it through a CLT. Direct use of Stein's method fails for submer sketches because usually they do not determine the sequence length exactly. In principle, however, Stein's method provides an approximate heuristic bound through estimation of the sequence length. In practice, the "estimate method" fails for many parameter values, because the sparsity of submers causes the bound from Stein's method to be too generous. Accordingly, our code relies on a "qualitative method", leaving simulations to estimate the accuracy of the confidence intervals from the CLT convergence.  

A. Dutta et al. (2022)<br/>
Parameterized syncmer schemes improve long-read mapping<br />
bioRxiv: 2022.01.10.475696

influenced us to generalize our methods for open and closed syncmers to their "parametrized syncmers".

Our code estimates sequence length and mutation probability per base by using submer counts in central limit theorems. Presently, the submers can be of 3 types: parametrized syncmers, minimizers, or minimally overlapping k-mers. Note, however, the context-dependency of minimizers obstructed the estimation of the corresponding mutation probabilities. The code also calculates the first-occurrence probabilities (the inter-submer distance distribution) for each submer type. A general formula converts the first-occurrence probabilities to the alpha-run probabilities of  

J. Shaw & Y.W. Yu (2021)<br />
Theory of local k-mer selection with applications to long-read alignment<br />
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab790/6432031<br />

**programs/** contains Python executables and ...make.py drivers, which display calls to the programs.<br />
**Output/** contains the output from the ...make.py drivers, so diff can verify the executables by comparing their output files ...log.<br />
**modules/** contains Python modules and classes, tested by a main() to demonstrate calls to the module subroutines.<br />

**programs/** 
Executables have an -h (help) option to explain their arguments. 
1. **distance-distribution...py** outputs first-occurrence probabilities. If the '-y' flag is set, it outputs alpha-test probabilities.
2. **length-from...py** inputs the submer count of a sequence and outputs a confidence interval for its length.
3. **theta-from...py** inputs the submer counts of a reference sequence and a mutated version and outputs a confidence interval for the mutation probability per letter.

**modules/** contains the files performing the computations. A brief summary of the most important of these files follow.

1. **jls_submer_clt_mgr.py** (manager file for submer CLTs) is the main programming interface. Please refer to its comments on the arguments and the return of its subroutines for more information on the CLTs.

2. For all types of submers, let the indicator Y<sub>i</sub> = 1 if the i-th k-mer is a submer, and 0 otherwise. The CLTs depend on the the autocovariance function cov[Y<sub>0</sub>,Y<sub>i</sub>]. The base class file **submer.py** calculates the autocovariance from the expected products E[Y<sub>0</sub>Y<sub>i</sub>] provided by the derived class files for each submer type (**jls_syncmer_parametrized.py**, **jls_minimizer.py**, and **jls_non_overlapping_pattern_prob.py**). 

**Some miscellaneous topics**

1. The &alpha;-test probabilities of Shaw & Yu

The derived class files (**jls_syncmer_parametrized.py**, **jls_syncmer_minimizer.py**, and **jls_non_overlapping_pattern_prob.py**) calculate first-occurrence probabilities (inter-submer distance distribution)

f<sub>i</sub> = Pr{ Y<sub>i</sub> = 1 and Y<sub>j</sub> = 0 (0 < j < i) | Y<sub>0</sub> = 1 }.

In contrast, the expected products E[Y<sub>0</sub>Y<sub>i</sub>] leave Y<sub>j</sub> for 0 < j < i unrestricted. 

**submer.py** interconverts the &alpha;-test probabilities Pr(f,&alpha;) of Yun & Shaw and first-passage probabilities f<sub>i</sub>, which determine each other according to relatively simple formulas. Yun & Shaw give general four-variable recursions for the &alpha;-test probabilities, but the recursions here are much faster.

2. Window probabilities for (w,k)-minimizers

Minimizers have a window guarantee: if two sequences share a (w+k-1)-mer, then they share a (w,k)-minimizer. Subsequences of length &alpha;+k-1 only have a window probability &pi;(&alpha;), where &pi;(&alpha;) = 1 for &alpha; &ge; w; and &pi;(&alpha;) &lt; 1 otherwise. **jls_minimizer_common_submer.py** calculates &pi;(&alpha;).



 
