# submer-central-limit-theorem

## The central limit theorems  (**CLT**s)

The code was inspired by the preprint

A. Blanca et al. (2021)<br/> 
The statistics of k-mers from a sequence undergoing a simple mutation process without spurious matches.<br/> 
DOI: 10.1101/2021.01.15.426881

which used Stein's method to develop CLTs for the k-mer sketch of sequences. Under reasonable assumptions, the k-mer sketch can determine sequence length exactly because the it lists all k-mers in the sequence. In contrast, submers sample the k-mers in a sequence. Thus, the submer count does not determine the sequence length exactly, but yields probabilisitic bounds on it through a CLT. Direct use of Stein's method fails for submer sketches because usually they do not determine the sequence length exactly. In principle, however, Stein's method provides an approximate heuristic bound through an "estimate method". In practice, the estimate method fails for many parameter values, because the sparsity of submers causes the bound from Stein's method to be too generous. Accordingly, the code relies on a "qualitative method", leaving simulations to estimate the speed of the CLT convergence.  

Primarily, the code estimates sequence length and mutation probability per base by using submer counts in central limit theorems. Presently, the submers can be of 4 types: closed syncmers, open syncmers, minimizers, or minimally overlapping k-mers. Note, however, the context-dependency of minimizers obstructed the estimation of the corresponding mutation probabilities. The code also calculates the first-passage probabilities (inter-submer distance distribution) for each submer type. A general formula converts the first-passage probabilities to the alpha-test probabilities of  

J. Shaw & Y.W. Yu (2021)<br />
Theory of local k-mer selection with applications to long-read alignment<br />
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab790/6432031<br />

## **The main directory** contains Python executables.

The executables have the -h (help) option to explain their arguments. The system calls in various *_make.py files display the arguments and the *.log files display the output of the executables.

1. **distance-distribution*.py** outputs first-passage probabilities. If the '-y' flag is set, it outputs alpha-test probabilities.

2. **length-from*.py** outputs a confidence interval for the sequence length from the submer count of a sequence.

3. **theta-from*.py** outputs a confidence interval for the mutation probability per letter from the submer count of a sequence and the submer counts common to a reference sequence and a mutated version.

## **modules/** contains the files performing the computations.

contains the computations. **jls_submer_clt_mgr.py** (manager file for submer CLTs) is the main programming interface. Please refer to its comments on the arguments and the return of its subroutines for the CLTs.

For all types of submers, the CLTs depend on the following autocovariance function. Let the indicator Y<sub>i</sub> = 1 if the i-th k-mer is a submer, and 0 otherwise. The base class file **submer.py** calculates cov[Y<sub>0</sub>,Y<sub>i</sub>] from the expected products E[Y<sub>0</sub>Y<sub>i</sub>] provided by the derived class files for each submer type (**jls_syncmer_closed.py**, **jls_syncmer_open.py**, **jls_minimizer.py**, and **jls_non_overlapping_pattern_prob.py**). 

## The &alpha;-test probabilities of Shaw & Yu

The derived class files (**jls_syncmer_closed.py**, **jls_syncmer_open.py**, **jls_syncmer_minimizer.py**, and **jls_non_overlapping_pattern_prob.py**) calculate the first-passage probabilities 

f<sub>i</sub> = Pr{ Y<sub>i</sub> = 1 and Y<sub>j</sub> = 0 (0 < j < i) | Y<sub>0</sub> = 1 }.

Note in contrast, the expected products E[Y<sub>0</sub>Y<sub>i</sub>] leave Y<sub>j</sub> for 0 < j < i unrestricted. 

The file **submer.py** interconverts the &alpha;-test probabilities Pr(f,&alpha;) of Yun & Shaw and first-passage probabilities f<sub>i</sub>, which always determine each other. Yun & Shaw give general four-variable recursions for the &alpha;-test probabilities, but the recursions here are much faster.

## Window probabilities for (w,k)-minimizers

Minimizers have a window guarantee: if two sequences share a (w+k-1)-mer, then they share a (w,k)-minimizer. Subsequences of length &alpha;+k-1 only have a window probability &pi;(&alpha;), where &pi;(&alpha;) = 1 for &alpha; &ge; w; and &pi;(&alpha;) &lt; 1 otherwise. **jls_minimizer_common_submer.py** calculates &pi;(&alpha;).



 
