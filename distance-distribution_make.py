#!/usr/bin/env python

from os import system

log = 'distance-distribution.log'
open(log, 'w').close() # The file {log} is empty.

kmer_length = 6
smer_length = 2
max_distance = 5
is_yu_and_shaw = ''

# Let u = k-s.

# density of closed syncmers = 2 / (u+1)
# Closed syncmers have a guarantee for a u-window of k-mers.
# Equation numbers refer to Word document on closed syncmers, 
# Tests first_passage_probabilities, f[i] by hand calculation.

# Eq (24) gives first-passage probability f[i].
# f = [0, 1/3, 4/21, 5/42, 5/14, 0] # The 4-window guarantee makes further computation pointless.
command = f'python distance-distribution-closed-syncmer.py -k {kmer_length} -s {smer_length} -m {max_distance} {is_yu_and_shaw}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

# Tests Yu & Shaw Pr vector, Pr[i] by hand calculation.
# Eq (28) gives Yu & Shaw Pr.
# Pr(f,alpha) = 2*alpha/(u+alpha)
# Pr(f,alpha) = [0, 2/5, 4/6, 6/7, 8/8, 1] 
is_yu_and_shaw = '-y'
command = f'python distance-distribution-closed-syncmer.py -k {kmer_length} -s {smer_length} -m {max_distance} {is_yu_and_shaw}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

kmer_length = 6
smer_length = 2
t_offset_of_smer = 4
max_distance = 23
is_yu_and_shaw = ''

# density of open syncmers = 1 / (u+1)
# Equation numbers refer to Word document on open syncmers, 
# Tests first_passage_probabilities, f[i] by hand calculation.

# Eq (24) gives first-passage probability f[i].
# f = [0, 1/3, 4/21, 5/42, 5/14, 0] 
command = f'python distance-distribution-open-syncmer.py -k {kmer_length} -s {smer_length} -t {t_offset_of_smer} -m {max_distance} {is_yu_and_shaw}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

# Tests Yu & Shaw Pr vector, Pr[i] by hand calculation.
# Eq (28) gives Yu & Shaw Pr.
# Pr(f,alpha) = 2*alpha/(u+alpha)
# Pr(f,alpha) = [0, 2/5, 4/6, 6/7, 8/8, 1] 
is_yu_and_shaw = '-y'
command = f'python distance-distribution-open-syncmer.py -k {kmer_length} -s {smer_length} -t {t_offset_of_smer} -m {max_distance} {is_yu_and_shaw}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")
