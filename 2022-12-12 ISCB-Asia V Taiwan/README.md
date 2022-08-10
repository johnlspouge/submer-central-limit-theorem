## Parametrized syncmers and their inter-syncmer distance distribution and run-hitting probabilities

The code in programs/distance-distribution-parametrized-syncmer.py reflects concepts in the preprint of

A. Dutta et al. (2022)<br/>
Parameterized syncmer schemes improve long-read mapping<br />
bioRxiv: 2022.01.10.475696

which generalizes open and closed syncmers to "parametrized syncmers". The calls below use the parametrization, e.g.,

python distance-distribution-parametrized-syncmer.py -k 6 -s 2 -t 0 4 -m 5<br />

The command finds the inter-syncmer distance distribution for (k,s,t)=(6,2,{0,4}) closed syncmers out to distance 5. (Closed syncmers actually have a window guarantee of 4.) Similarly,

python distance-distribution-parametrized-syncmer.py -k 6 -s 2 -t 1 -m 5 -y<br />

finds the inter-syncmer distance distribution for (k,s,t)=(6,2,{1}) open syncmers out to distance 23. Open syncmers lack a window guarantee. 

Parametrized syncmers also include an option -e for downsampling, e.g., adding an option' -e 0.1' indicates a rejection probability eps=0.1 for downsampling. The examples omit the option and default to eps=0.0, the corresponding values with no syncmer downsampling.<br />

A general formula given in the conference paper converts the inter-syncmer distance distribution to the run-hitting probabilities of  

J. Shaw & Y.W. Yu (2021)<br />
Theory of local k-mer selection with applications to long-read alignment<br />
https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btab790/6432031<br />

**programs/** contains the Python executable<br /> 
**distance-distribution.py** The executable has an -h (help) option to explain its arguments.<br /> 
which outputs the inter-syncmer distance distribution. If the '-y' flag is set, it outputs run-hitting probabilities.distance-distribution-parametrized-syncmer.py. The driver 
**distance-distribution_make.py** contains calls to the executable.<br />
**programs/Output/** contains the output from the calls, so<br />
diff distance-distribution.log programs/Output/<br />
can verify distance-distribution.py by comparing the output file distance-distribution.log to a copy in programs/Output/.<br />

**modules/** contains the files performing the computations. The python driver test_make.py calls unit tests in each file. Each call exits without complaint if the result is consistent with previous computations.
