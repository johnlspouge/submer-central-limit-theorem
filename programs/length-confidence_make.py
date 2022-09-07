#!/usr/bin/env python

from os import system

log = 'length-confidence.log'
open(log, 'w').close() # The file {log} is empty.

counts = [100, 10000, 1000000]
flags = [''] # The ' -e' flag never produced anything but [None, None] for me.
submer0params = (('parametrized-syncmer',' -k 10 -s 3 -t 0 7'), # density 2/(10-3+1) = 0.25
                 ('parametrized-syncmer',' -k 10 -s 3 -t 1'),  # density 1/(10-3+1) = 0.125
                 ('minimizer',' -k 10 -w 3'), 
                 ('non-overlapping-pattern',' -k 6 -p "a,cg" -f "0.15,0.2,0.5,0.15"')) # density 0.25

for flag in flags:
    for submer0param in submer0params:
        key,value = submer0param
        for count in counts:
            command = f'python length-confidence-from-{key}-count.py {value} -n {count} -c 0.95 {flag}'
            with open(log, "a") as fh:
                fh.write(command + "\n")
            command = f'{command} >> {log}'
            system( f'{command}' )
            with open(log, "a") as fh:
                fh.write("\n")
