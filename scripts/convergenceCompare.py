#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

pathToSimFolder1 = argv[1]
pathToSimFolder2 = argv[2]


with h5py.File(join(pathToSimFolder1, "data.hdf5"), "r") as dataFile:
    optEnergy1 = np.array(dataFile["/optEnergy"][:])
with h5py.File(join(pathToSimFolder2, "data.hdf5"), "r") as dataFile:
    optEnergy2 = np.array(dataFile["/optEnergy"][:])


minIterations = min(optEnergy1.shape[0], optEnergy2.shape[0])


fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))


ax.plot(
    np.maximum.accumulate(optEnergy1[:minIterations]),
    "-",
    ms=1,
    color="darkred",
    label="genetic",
)
ax.plot(
    np.maximum.accumulate(optEnergy2[:minIterations]),
    "-",
    ms=1,
    color="darkblue",
    label="MC",
)


# ax.set_xlim(-0.15,0.65)
# ax.set_ylim(-1.1,1.1)

ax.set_xlabel(r"iteration")
ax.set_ylabel(r"$\mathcal F$")


ax.legend()


plt.savefig(
    join(pathToSimFolder1, "convergenceCompare.png"), bbox_inches="tight", dpi=300
)
# plt.savefig(join(pathToSimFolder2,"convergenceCompare.png"),bbox_inches="tight",dpi=300)
plt.show()
plt.close(fig)
