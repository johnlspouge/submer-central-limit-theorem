# submer-central-limit-theorem

## Parametrized syncmers

The code reflects concepts in the preprint of

A. Dutta et al. (2022)<br/>
Parameterized syncmer schemes improve long-read mapping<br />
bioRxiv: 2022.01.10.475696

which generalized open and closed syncmers to "parametrized syncmers". The calls below use parameters for syncmers, e.g.,

python distance-distribution-parametrized-syncmer.py -k 6 -s 2 -t 0 4 -m 5<br />

finds the inter-syncmer distance distribution for (k,s,t)=(6,2,{0,4}) (closed) syncmers out to distance 5, where closed syncmers have a window guarantee of 4. Similarly,

python distance-distribution-parametrized-syncmer.py -k 6 -s 2 -t 1 -m 5 -y<br />

finds the inter-syncmer distance distribution for open (k,s,t)=(6,2,{1}) syncmers out to distance 23. Open syncmers lack a window guarantee. 

where closed (k,s)=(10,3) syncmers are now parametrized (k,s, ts)=(10, 3,'0 7') syncmers,<br />
and open (k,s,t)=(10,3,1) syncmers are now parametrized (k,s, ts)=(10, 3,'1') syncmers.

Parametrized syncmers now include an option -e for downsampling, e.g., ' -e 0.1' indicates a downsampling probability eps=0.1.<br />
Omit the option to default to eps=0.0, the corresponding values with no syncmer downsampling.

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

**submer.py** interconverts the &alpha;-test probabilities Pr(f,&alpha;) of Shaw & Yu (2021) and first-passage probabilities f<sub>i</sub>, which determine each other according to relatively simple formulas. Shaw & Yu (2021) give general four-variable recursions for the &alpha;-test probabilities, and Dutta et al (2022) other recursions, but the recursions given here are much faster.

2. Window probabilities for (w,k)-minimizers

Minimizers have a window guarantee: if two sequences share a (w+k-1)-mer, then they share a (w,k)-minimizer. Subsequences of length &alpha;+k-1 only have a window probability &pi;(&alpha;), where &pi;(&alpha;) = 1 for &alpha; &ge; w; and &pi;(&alpha;) &lt; 1 otherwise. **jls_minimizer_common_submer.py** calculates &pi;(&alpha;).



 
