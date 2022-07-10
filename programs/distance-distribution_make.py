#!/usr/bin/env python

from os import system

log = 'distance-distribution.log'
open(log, 'w').close() # The file {log} is empty.

# Let u = k-s.

# density of closed syncmers = 2 / (u+1)
# Closed syncmers have a guarantee for a u-window of k-mers.
# Equation numbers refer to Word document on closed syncmers. 

# Run 0
# Tests first_occurrence_probabilities, f[i] by hand calculation.
# Eq (24) gives first-occurrence probability f[i].
# f = [0, 1/3, 4/21, 5/42, 5/14, 0] # The 4-window guarantee makes further computation pointless.

# Run 1
# Tests Yu & Shaw Pr vector, Pr[i] by hand calculation.
# Eq (28) gives Yu & Shaw Pr.
# Pr(f,alpha) = 2*alpha/(u+alpha)
# Pr(f,alpha) = [0, 2/5, 4/6, 6/7, 8/8, 1] 

executable = 'distance-distribution-parametrized-syncmer.py'
kmer_length = 6
smer_length = 2
t_offsets_of_smer = '0 4'
max_distance = 5
is_test_probabilities_runs = ['', '-y']

for is_test_probabilities in is_test_probabilities_runs:
    command = f'python {executable} -k {kmer_length} -s {smer_length} -t {t_offsets_of_smer} -m {max_distance} {is_test_probabilities}'
    with open(log, "a") as fh:
        fh.write(command + "\n")
    command = f'{command} >> {log}'
    system( f'{command}' )
    with open(log, "a") as fh:
        fh.write("\n")

# density of open syncmers = 1 / (u+1)
# Equation numbers refer to Word document on open syncmers. 

# Run 0
# Tests first_occurrence_probabilities, f[i] by hand calculation.
# Eq (24) gives first-occurrence probability f[i].

# Run 1
# Tests Yu & Shaw Pr vector, Pr[i] by hand calculation.
# Eq (28) gives Yu & Shaw Pr.

executable = 'distance-distribution-parametrized-syncmer.py'
kmer_length = 6
smer_length = 2
t_offsets_of_smer = '4'
max_distance = 23
is_test_probabilities_runs = ['', '-y']

for is_test_probabilities in is_test_probabilities_runs:
    command = f'python {executable} -k {kmer_length} -s {smer_length} -t {t_offsets_of_smer} -m {max_distance} {is_test_probabilities}'
    with open(log, "a") as fh:
        fh.write(command + "\n")
    command = f'{command} >> {log}'
    system( f'{command}' )
    with open(log, "a") as fh:
        fh.write("\n")

# density of minimally overlapping k-mers = P(H)
# Equation numbers refer to Word document on minimally overlapping k-mers. 

# Run 0
# Tests first_occurrence_probabilities, f[i] by hand calculation.
# Eq (24) gives first-occurrence probability f[i].

# Run 1
# Tests Yu & Shaw Pr vector, Pr[i] by hand calculation.
# Eq (28) gives Yu & Shaw Pr.

executable = 'distance-distribution-non-overlapping-pattern.py'
kmer_length = 6
patterns = 'a,cg'
frequencies = '0.15,0.2,0.5,0.15'
max_distance = 23
is_test_probabilities_runs = ['', '-y']

for is_test_probabilities in is_test_probabilities_runs:
    command = f'python {executable} -k {kmer_length} -p {patterns} -f {frequencies} -m {max_distance} {is_test_probabilities}'
    with open(log, "a") as fh:
        fh.write(command + "\n")
    command = f'{command} >> {log}'
    system( f'{command}' )
    with open(log, "a") as fh:
        fh.write("\n")