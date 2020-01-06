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

with h5py.File(join(pathToSimFolder,"data.hdf5"),"r") as dataFile:
    currents=np.array(dataFile["/outputCurrent"][:])
    sigma=np.array(dataFile["/outputCurrentUncert"][:])
    fitness=np.array(dataFile["/fitness"][:])
    voltages=np.array(dataFile["/voltages"][:])

voltageScanPointsNumber = int((parameters["voltageScanMax"]-parameters["voltageScanMin"])/parameters["voltageScanResolution"]+1)
currents.resize((currents.shape[0],voltageScanPointsNumber,voltageScanPointsNumber))
sigma   .resize((sigma   .shape[0],voltageScanPointsNumber,voltageScanPointsNumber))

print(currents.shape)

best=np.argmax(fitness,axis=0)
print("best at",best,"fitness: ",fitness[best]," +-", sigma[best])
print(voltages[best])


fig, ax = plt.subplots(1,1)

ax.errorbar(np.arange(4),currents[best].flatten(),yerr=sigma[best].flatten(),marker="x",color="k")

plt.show()



# fig, ax = plt.subplots(1,1)

# N=10
# for i in range(N):
#     ax.plot(np.arange(4),currents[i].flatten(),"x",color=color(i,N))
# plt.show()



# fig, ax = plt.subplots(1,1)

# im=ax.imshow(currents[0],cmap="jet")

# plt.colorbar(im)

# plt.show()

