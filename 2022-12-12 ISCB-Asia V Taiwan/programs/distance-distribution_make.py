#!/usr/bin/env python

from os import system

log = 'distance-distribution.log'
open(log, 'w').close() # The file {log} is empty.

# Let u = k-s.

# density of closed syncmers = 2 / (u+1)
# Closed syncmers have a guarantee for a u-window of k-mers.

# Run ''
# Tests probability distribution of inter-submer distances.

# Run '-y'
# Tests Yu & Shaw Pr vector.

executable = 'distance-distribution-parametrized-syncmer.py'
kmer_length = 6
smer_length = 2
t_offsets_of_smer = '0 4' # 0-offset, not 1-offset
max_distance = 5 # maximum distance to be tests
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

# Run ''
# Tests probability distribution of inter-submer distances.

# Run '-y'
# Tests Yu & Shaw Pr vector.

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
