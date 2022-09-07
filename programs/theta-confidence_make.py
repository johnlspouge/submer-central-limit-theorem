#!/usr/bin/env python

from os import system

log = 'theta-confidence.log'
open(log, 'w').close() # The file {log} is empty.

counts = [100, 1000, 10000]
submer0params = (('parametrized-syncmer',' -k 10 -s 3 -t 0 7'), # closed syncmer density 0.25
                 ('parametrized-syncmer',' -k 10 -s 3 -t 1'), # open syncmer density 0.125
                 ('non-overlapping-pattern',' -k 6 -p "a,cg" -f "0.15,0.2,0.5,0.15"')) # density 0.25
# without length estimate
#for submer0param in submer0params:
for submer0param in submer0params:
    submer,param = submer0param
    for count in counts:
        unmutated_count = int( 0.9 * count + 0.5 )
        command = f'python theta-confidence-from-{submer}-count.py {param} -n {count} -u {unmutated_count}  -c 0.95'
        with open(log, "a") as fh:
            fh.write(command + "\n")
        command = f'{command} >> {log}'
        system( f'{command}' )
        with open(log, "a") as fh:
            fh.write("\n")

# with length estimate
submer,param = submer0params[0]
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 4
command = f'python theta-confidence-from-{submer}-count.py {param} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

submer,param = submer0params[1]
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 8
command = f'python theta-confidence-from-{submer}-count.py {param} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

submer,param = submer0params[2]
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 4
command = f'python theta-confidence-from-{submer}-count.py {param} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")
