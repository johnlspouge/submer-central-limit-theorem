#!/usr/bin/env python

from os import system

log = 'length-hypothesis.log'
open(log, 'w').close() # The file {log} is empty.

submer0params = (
    ('parametrized-syncmer',' -k 10 -s 3 -t 0 7',0.25), 
    # density 2/(10-3+1) = 0.25
    ('parametrized-syncmer',' -k 10 -s 3 -t 1',0.125),  
    # density 1/(10-3+1) = 0.125
)

counts = [1000,10000,100000]
for count in counts:
    for submer0param in submer0params:
        submer,options,density = submer0param
        for i in range(2):
            n = int(count*(density+0.005*(-1)**(i+1)))
            command = f'python length-hypothesis-from-{submer}-count.py {options} -n {n} -0 {count}'
            with open(log, "a") as fh:
                fh.write(command + "\n")
            command = f'{command} >> {log}'
            system( f'{command}' )
            with open(log, "a") as fh:
                fh.write("\n")
