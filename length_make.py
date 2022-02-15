#!/usr/bin/env python

from os import system

log = 'length.log'
open(log, 'w').close() # The file {log} is empty.

counts = [100, 10000, 1000000]
flags = [''] # The ' -e' flag never produced anything but [None, None] for me.
submer2params = {'closed-syncmer':' -k 10 -s 3', 
                 'open-syncmer':' -k 10 -s 3 -t 1', 
                 'minimizer':' -k 10 -w 3'}

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
