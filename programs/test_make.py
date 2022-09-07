#!/usr/bin/env python

from os import system, listdir

dir = './'

# Runs the files *_make.py. 
for file in listdir(dir):
    if file.endswith("_make.py") and file != __file__:
        print(file)
        system( f'python {dir}{file}')
