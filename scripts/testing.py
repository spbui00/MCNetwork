#!/usr/bin/python3

from os.path import join
from re import sub
import numpy as np

def calcMemory(electrodes, acceptors):
    states=2**acceptors
    mem=((acceptors+electrodes)**2+1)*states*8
    print(f"need max mem of {mem/10**9} GB for {states} states")
    
calcMemory(8,30)