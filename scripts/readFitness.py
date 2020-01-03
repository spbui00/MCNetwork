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

currents=np.array(dataFile["/outputCurrent"][:])

voltageScanPointsNumber = int((parameters["voltageScanMax"]-parameters["voltageScanMin"])/parameters["voltageScanResolution"]+1)
currents.resize((voltageScanPointsNumber,voltageScanPointsNumber))

print(currents)

fig, ax = plt.subplots(1,1)


ax.imshow(currents,cmap="jet")


# plt.colorbar()

plt.show()
