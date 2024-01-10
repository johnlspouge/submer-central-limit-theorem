#!/usr/bin/env python

from os import system

log = 'theta-hypothesis.log'
open(log, 'w').close() # The file {log} is empty.

submer0params = (
    ('parametrized-syncmer',' -k 10 -s 3 -t 0 7',0.25,10), 
    # density 2/(10-3+1) = 0.25
    ('parametrized-syncmer',' -k 10 -s 3 -t 1',0.125,10),  
    # density 1/(10-3+1) = 0.125
)

counts = [100,1000,10000]

for c,count in enumerate(counts):
    unmutated_count = int( 0.9 * count + 0.5 )
    for submer0param in submer0params:
        submer,options,density,k = submer0param
        for i in range(2):
            for j in range(2):
                if j == 0 and c == 0:
                    estimate = ''
                elif j == 0:
                    continue
                else:
                    cod = count/density
                    estimate = f' -l {cod}'
                theta_mean = 1.0 - (unmutated_count/count)**(1.0/k)
                theta0 = theta_mean+0.0005*(-1)**(i+1)
                command = f'python theta-hypothesis-from-{submer}-count.py {options} -n {count} -u {unmutated_count} -0 {theta0} {estimate}'
                with open(log, "a") as fh:
                    fh.write(command + "\n")
                command = f'{command} >> {log}'
                system( f'{command}' )
                with open(log, "a") as fh:
                    fh.write("\n")
