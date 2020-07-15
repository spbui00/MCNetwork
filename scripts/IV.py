#!/usr/bin/python3
from tools import readParameters
from sys import argv

import h5py
import matplotlib.pylab as plt

import numpy as np
import os

if len(argv) > 1:
    pathToSimFolder = argv[1]
else:
    pathToSimFolder = "data/"

parameters, electrodes = readParameters(pathToSimFolder)
electrodeNumber = len(electrodes)

data = h5py.File(os.path.join(pathToSimFolder, "data.hdf5"), "r")
voltages = np.linspace(
    parameters["voltageScanMin"],
    parameters["voltageScanMax"],
    int(parameters["voltageScanPoints"]),
)
currents = data["outputCurrent"][0]

plt.plot(voltages, currents)
plt.show()
