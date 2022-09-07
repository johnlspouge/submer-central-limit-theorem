#!/usr/bin/env python

from os import system

log = 'island-statistics.log'
open(log, 'w').close() # The file {log} is empty.

dir = './'

# Eq (27) gives "islands:genome_fraction_of_starts" mu_y * sum f[k+alpha].
# Eq (28) gives "islands:genome_fraction_of_starts" mu_y * sum alpha * f[k+alpha].
# Eq (29) gives "islands:genome_fraction_of_starts" mu_y * sum alpha**2 * f[k+alpha].

# Let u = k-s.

# density of closed syncmers = 2 / (u+1)
# Closed syncmers have a guarantee for a u-window of k-mers.

executable = 'island-statistics-parametrized-syncmer.py'
kmer_length = 6
smer_length = 2
t_offsets_of_smer = '0 4'

command = f'python {dir}{executable} -k {kmer_length} -s {smer_length} -t {t_offsets_of_smer}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

# density of open syncmers = 1 / (u+1)

executable = 'island-statistics-parametrized-syncmer.py'
kmer_length = 6
smer_length = 2
t_offsets_of_smer = '4'

command = f'python {dir}{executable} -k {kmer_length} -s {smer_length} -t {t_offsets_of_smer}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

# density of minimally overlapping k-mers = P(H)

executable = 'island-statistics-non-overlapping-pattern.py'
kmer_length = 6
patterns = 'a,cg'
frequencies = '0.15,0.2,0.5,0.15'

command = f'python {dir}{executable} -k {kmer_length} -p {patterns} -f {frequencies}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")
