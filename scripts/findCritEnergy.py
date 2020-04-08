#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

from time import sleep

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"




runs = 3

fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

for r in range(runs):
    print(r+1)
    fileOpenTries = 0
    while fileOpenTries < 50:
        fileOpenTries += 1
        try:
            with h5py.File(join(pathToSimFolder,f"run_{r+1}/data.hdf5"),"r") as dataFile:
                optEnergy=np.array(dataFile["/optEnergy"][:,0])
            break
        except OSError as e:
            if "No such file" in repr(e) :
                raise e
            else:
                print(f"could not open file. try number {fileOpenTries}")
                sleep(1)
                

    accumulated = np.maximum.accumulate(optEnergy)

    distance = 50

    dataX = []
    dataY = []

    currentValue = 0
    shiftedValue = 0

    N = accumulated.shape[0]
    for i in range(N):
        if accumulated[i] != currentValue or accumulated[min(i+distance, N-1)] != shiftedValue:
            currentValue = accumulated[i]
            shiftedValue = accumulated[min(i+distance, N-1)]
            
            dataX.append(currentValue)
            dataY.append(shiftedValue)
            
            

    fig2, ax2=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


    ax2.plot(dataX, dataY,"k-")


    ax2.set_xlim(0.45,1)
    ax2.set_ylim(0.45,1)

    ax2.set_xlabel(r"$\mathcal{F}_{i}$")
    ax2.set_ylabel(rf"$\mathcal{{F}}_{{i+{distance} }}$")


    plt.savefig(join(pathToSimFolder,f"run_{r+1}/critFitness_{r+1}.png"),bbox_inches="tight",dpi=300)
    # plt.show()
    plt.close()
    fig=None
    
    
    ax.plot(dataX, dataY,"k-",color = color(r,runs))
    




ax.set_xlim(0.45,1)
ax.set_ylim(0.45,1)

ax.set_xlabel(r"$\mathcal{F}_{i}$")
ax.set_ylabel(rf"$\mathcal{{F}}_{{i+{distance} }}$")


plt.savefig(join(pathToSimFolder,"critFitness.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close()
fig=None