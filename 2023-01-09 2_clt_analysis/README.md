# syncmer-central-limit-theorem

## The central limit theorems  (CLTs)

The contents of this directory implement confidence intervals and hypothesis tests for parametrized syncmers (with possible downsampling).

### 1. Contents of programs/

The 'test_make.py' file should run and exit quietly within about 5 minutes if the Python executables are functioning as expected. Besides 'test_make.py', programs/ contains 12 other files relevant to CLTs for parametrized syncmers. The 12 filenames are a concatenation of 4 [PREFIX]-es and 3 [POSTFIX]-es.

The 4 [PREFIX]-es are:

    A.  'length-confidence'
    B.  'length-hypothesis'
    C.  'theta-confidence'
    D.  'theta-hypothesis'

The 3 [POSTFIX]-es are:

    A. '-from-parametrized-syncmer-count.py':
        The executable implements the confidence interval or hypothesis test above.
    B. '_make.py':
        The makefile tests the implementations of the confidence interval or hypothesis test.
    C. '.log':
        The file is output from the corresponding '_make.py' makefile. 
        Pairs of successive lines display: 
            i. an example of the program calling for a confidence interval or hypothesis test 
            ii. the corresponding program output.

The calls to executables in A. above have forms like
* 'python length-confidence-from-parametrized-syncmer-count.py [ARGUMENTS]'

The argparse module lists and explains the [ARGUMENTS] with a -h (help) call of the form 
* 'python length-confidence-from-parametrized-syncmer-count.py -h'

The following provide examples of calls to the executables in A. above, followed by their output.

* 'python length-confidence-from-parametrized-syncmer-count.py  -k 10 -s 3 -t 0 7 -n 100 -c 0.95'
* sig	method	LengthLow	LengthHigh<br/>
  0.95	qualitative	367.3508315817253	435.4890933863735

*Explanation*: Closed (k,s)=(10,3) syncmers are parametrized as (k, s, ts)=(10, 3, '0 7') syncmers, placing the s-minimizer at the beginning or end of the k-mer. **All calls to our parametrized syncmer code use offset from array index 0 within the k-mer**, thus ts = '0 7' are the offsets of the initial and terminal s-mers in a closed syncmer. The total number of syncmers is -n 100; and the confidence level is -c 0.95 (i.e., 95%). The 95% confidence interval for the sequence length is (367.3508315817253,	435.4890933863735). 

* python length-hypothesis-from-parametrized-syncmer-count.py  -k 10 -s 3 -t 1 -n 120 -0 1000
* p_1-sided_left	p_1-sided_right	p_2-sided
  0.18002805028895036	0.8199719497110496	0.3600561005779007

*Explanation*: An open (k, s, t)=(10, 3, 1) places the s-minimizer at the offset after the beginning of the k-mer. **All calls to our parametrized syncmer code use offset from array index 0 within the k-mer**, thus ts = '1' is the offset after the initial s-mers in the open syncmer. The total number of syncmers is -n 120; and the hypothesized sequence length is -0 1000. The (1-sided_left p-value) probability that the length is as small as 1000 or less is 0.18002805028895036; (1-sided_right p-value) probability that the length is as large as 1000 or greater is 0.8199719497110496; and the (2-sided p-value) probability that the length is at least as extreme as 1000 is 0.3600561005779007 = 2.0*min(0.18002805028895036,	0.8199719497110496). 

* 'python theta-confidence-from-parametrized-syncmer-count.py  -k 10 -s 3 -t 0 7 -n 100 -u 90  -c 0.95'
* sig	length	ThetaLow	ThetaHigh<br/>
  0.95	estimated	0.0038246028859374994	0.026768396719205736

*Explanation*: Closed (k,s)=(10,3) syncmers are parametrized as (k, s, ts)=(10, 3, '0 7') syncmers, placing the s-minimizer at the beginning or end of the k-mer. **All calls to our parametrized syncmer code use offset from array index 0 within the k-mer**, thus ts = '0 7' are the offsets of the initial and terminal s-mers in the closed syncmer. The total number of syncmers in sequence *A* is -n 100; the total number of unmutated syncmers between sequences *A* and *B* is -u 90; and the confidence level is -c 0.95 (i.e., 95%). The reference sequence length *A* was estimated from the total syncmer count of 100, yielding a 95% confidence interval (0.0038246028859374994,	0.026768396719205736) for the mutation probability &theta;. If the actual reference sequence length *A* is known, it can be added as a parameter -l [LENGTH].

* python theta-hypothesis-from-parametrized-syncmer-count.py  -k 10 -s 3 -t 0 7 -n 100 -u 90 -0 0.010980741793785553  -l 400.0
* length	p_1-sided_left	p_1-sided_right	p_2-sided
  actual	0.5354548492893464	0.4645451507106536	0.9290903014213072

*Explanation*: An closed (k, s, t)=(10, 3, '0 7') places the s-minimizer at the beginning or end of the k-mer. **All calls to our parametrized syncmer code use offset from array index 0 within the k-mer**, thus ts = '0 7' are the offsets of the initial and terminal s-mers in the closed syncmer. The total number of syncmers is -n 100; the length of the actual reference sequence is -l 400.0 (if omitted, the length is estimated from the total number of syncmers); and the hypothesized mutation probability &theta; is -0 0.010980741793785553. The (1-sided_left p-value) probability that &theta is as small as 0.010980741793785553 or less is 0.5354548492893464; (1-sided_right p-value) probability that &theta is as large as 0.010980741793785553 or greater is 0.4645451507106536; and the (2-sided p-value) probability that &theta is at least as extreme as 0.010980741793785553 is 0.9290903014213072 = 2.0*min(0.5354548492893464,	0.4645451507106536). 

### 2. Contents of modules/
The 'test_make.py' file should run and exit quietly within about 1 minute if the Python executables are functioning as expected. The calls in programs/ do not require knowledge of the contents of modules/. The directory modules/ contains files for implementing derived classes and utilities whose the elementary operations model downsampled parametrized syncmers. The hierarchy of derived classes in the files is Submer => Syncmer => Parametrized Syncmer.  


 
