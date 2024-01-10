#!/usr/bin/env python

from os import system, listdir
from os.path import basename

dir = './'

print(basename(__file__))

# Runs the files *_make.py. 
for file in sorted(listdir(dir)):
    if file.endswith("_make.py") and file != basename(__file__):
        print(file)
        system( f'python {dir}{file}')
