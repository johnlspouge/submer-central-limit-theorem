#!/usr/bin/env python

from os import system

log = 'theta.log'
open(log, 'w').close() # The file {log} is empty.

counts = [100, 1000, 10000]
submer2params = {'closed-syncmer':' -k 10 -s 3', # density 0.25
                 'open-syncmer':' -k 10 -s 3 -t 1', # density 0.125
                 'non-overlapping-pattern':' -k 6 -p "a,cg" -f "0.15,0.2,0.5,0.15"'} # density 0.25
# without length estimate
for key,value in submer2params.items():
    for count in counts:
        unmutated_count = int( 0.9 * count + 0.5 )
        command = f'python theta-from-{key}-count.py {value} -n {count} -u {unmutated_count} -c 0.95'
        with open(log, "a") as fh:
            fh.write(command + "\n")
        command = f'{command} >> {log}'
        system( f'{command}' )
        with open(log, "a") as fh:
            fh.write("\n")

# with length estimate
value = ' -k 10 -s 3'
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 4
command = f'python theta-from-closed-syncmer-count.py {value} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

value = ' -k 10 -s 3 -t 1'
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 8
command = f'python theta-from-open-syncmer-count.py {value} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")

value = ' -k 6 -p "a,cg" -f "0.15,0.2,0.5,0.15"'
count = 100
unmutated_count = int( 0.9 * count + 0.5 )
length_estimate = 100 * 4
command = f'python theta-from-non-overlapping-pattern-count.py {value} -n {count} -u {unmutated_count} -c 0.95 -l {length_estimate}'
with open(log, "a") as fh:
    fh.write(command + "\n")
command = f'{command} >> {log}'
system( f'{command}' )
with open(log, "a") as fh:
    fh.write("\n")
