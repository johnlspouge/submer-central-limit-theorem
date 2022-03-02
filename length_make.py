#!/usr/bin/env python

from os import system

log = 'length.log'
open(log, 'w').close() # The file {log} is empty.

counts = [100, 10000, 1000000]
flags = [''] # The ' -e' flag never produced anything but [None, None] for me.
submer2params = {'closed-syncmer':' -k 10 -s 3', # density 2/(10-3+1) = 0.25
                 'open-syncmer':' -k 10 -s 3 -t 1',  # density 1/(10-3+1) = 0.125
                 'minimizer':' -k 10 -w 3', 
                 'non-overlapping-pattern':' -k 6 -p "a,cg" -f "0.15,0.2,0.5,0.15"'} # density 0.25

for flag in flags:
    for key,value in submer2params.items():
        for count in counts:
            command = f'python length-from-{key}-count.py {value} -n {count} -c 0.95 {flag}'
            with open(log, "a") as fh:
                fh.write(command + "\n")
            command = f'{command} >> {log}'
            system( f'{command}' )
            with open(log, "a") as fh:
                fh.write("\n")
