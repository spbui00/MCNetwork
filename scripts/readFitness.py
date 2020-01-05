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
fitness=np.array(dataFile["/fitness"][:])

voltageScanPointsNumber = int((parameters["voltageScanMax"]-parameters["voltageScanMin"])/parameters["voltageScanResolution"]+1)
currents.resize((currents.shape[0],voltageScanPointsNumber,voltageScanPointsNumber))

print(currents.shape)
fig, ax = plt.subplots(1,1)

# for i in range(20):
    # print(i)
    # # ax.plot(np.arange(i*5,(i+1)*5),fitness[i*5:(i+1)*5],"x",color=color(i,20))
    # ax.errorbar([i*5+2.5],[np.mean(fitness[i*5:(i+1)*5])],yerr=np.std(fitness[i*5:(i+1)*5]),marker=".",color=color(i,20))

im=ax.imshow(currents[0],cmap="jet")

plt.colorbar(im)

plt.show()
