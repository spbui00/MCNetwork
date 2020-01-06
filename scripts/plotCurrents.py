#!/usr/bin/python3
from tools import readParameters, color

from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)

dataFile=h5py.File(join(pathToSimFolder,"data.hdf5"),"r")

currents=np.array(dataFile["/absCurrents"][:])


print(currents)

# fig, ax = plt.subplots(1,1)

# for i in range(currents.shape[1]):
#     ax.plot(np.arange(currents.shape[0]),currents[:,i],color=color(i,currents.shape[1]),label=f"{i}")


# ax.legend()

# plt.show()
