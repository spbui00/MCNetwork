#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

from time import sleep

if len(argv) > 1:
    pathToSimFolder = argv[1]
else:
    pathToSimFolder = "../data/"

parameters, electrodes = readParameters(pathToSimFolder)

fileOpenTries = 0
while fileOpenTries < 50:
    fileOpenTries += 1
    try:
        with h5py.File(join(pathToSimFolder, "data.hdf5"), "r") as dataFile:
            currents = np.array(dataFile["/outputCurrent"][:])
            sigma = np.array(dataFile["/outputCurrentUncert"][:])
            fitness = np.array(dataFile["/fitness"][:])
            sigmaFitness = np.array(dataFile["/fitnessUncert"][:])
            voltages = np.array(dataFile["/voltages"][:])
            optEnergy = np.array(dataFile["/optEnergy"][:])
        break
    except OSError as e:
        if "No such file" in repr(e):
            raise e
        else:
            print(f"could not open file. try number {fileOpenTries}")
            sleep(1)


voltageScanPoints = int(parameters["voltageScanPoints"])

resolution = (parameters["voltageScanMax"] - parameters["voltageScanMin"]) / (
    voltageScanPoints - 1
)


best = np.argmax(optEnergy, axis=0)
# best=np.argmax(fitness,axis=0)
# best=[0]
print(
    "best at",
    best,
    "optEnergy: ",
    optEnergy[best],
    "fitness: ",
    fitness[best],
    " +-",
    sigmaFitness[best],
)
print(voltages[best])

# print(currents[best])


if voltageScanPoints == 2:
    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.errorbar(
        np.arange(4),
        currents[best].flatten(),
        yerr=sigma[best].flatten(),
        marker="x",
        color="k",
        linestyle="",
    )

    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels([r"[0,0]", r"[1,0]", r"[0,1]", r"[1,1]"])

    ax.set_xlabel("gate input")
    ax.set_ylabel(r"$I^{\textrm{out}} [e \nu_0]$")

    plt.savefig(join(pathToSimFolder, "fitness_1D.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)


falseIndex = round((0 - parameters["voltageScanMin"]) / resolution)
trueIndex = round((0.5 - parameters["voltageScanMin"]) / resolution)


maxTrue = -np.inf
minTrue = np.inf

currents2x2 = np.zeros((2, 2))
currentUncert2x2 = np.zeros((2, 2))

for b1 in [0, 1]:
    for b2 in [0, 1]:
        if b1:
            i1 = trueIndex
        else:
            i1 = falseIndex
        if b2:
            i2 = trueIndex
        else:
            i2 = falseIndex

        currents2x2[b1, b2] = currents[best[0], i1, i2]
        currentUncert2x2[b1, b2] = sigma[best[0], i1, i2]

        # print(b1,b2,logical(b1,b2,parameters["gate"]))
        if logical(b1, b2, parameters["gate"]):
            if maxTrue < currents[best[0], i1, i2]:
                maxTrue = currents[best[0], i1, i2]
        else:
            if minTrue > currents[best[0], i1, i2]:
                minTrue = currents[best[0], i1, i2]

print(currents2x2)

fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

im = ax.imshow(
    (currents[best[0], ::-1, :] - minTrue) / (maxTrue - minTrue),
    cmap="Spectral",
    extent=[
        parameters["voltageScanMin"] - resolution / 2,
        parameters["voltageScanMax"] + resolution / 2,
        parameters["voltageScanMin"] - resolution / 2,
        parameters["voltageScanMax"] + resolution / 2,
    ],
)

ax.scatter([0, 0, 0.5, 0.5], [0, 0.5, 0, 0.5], c="r", marker="x")

# ax.set_xlim(-0.15,0.65)
# ax.set_ylim(-0.15,0.65)

ax.set_xlabel(r"$V^\textrm{in}_1$ [V]")
ax.set_ylabel(r"$V^\textrm{in}_2$ [V]")


plt.colorbar(im)

plt.savefig(join(pathToSimFolder, "fitness_2D.png"), bbox_inches="tight", dpi=300)
# plt.show()
plt.close(fig)


fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

im = ax.imshow(
    (currents[best[0], ::-1, :] - minTrue) / (maxTrue - minTrue),
    cmap="Spectral",
    vmin=0,
    vmax=1,
    extent=[
        parameters["voltageScanMin"] - resolution / 2,
        parameters["voltageScanMax"] + resolution / 2,
        parameters["voltageScanMin"] - resolution / 2,
        parameters["voltageScanMax"] + resolution / 2,
    ],
)

ax.scatter([0, 0, 0.5, 0.5], [0, 0.5, 0, 0.5], c="r", marker="x")

# ax.set_xlim(-0.15,0.65)
# ax.set_ylim(-0.15,0.65)

ax.set_xlabel(r"$V^\textrm{in}_1$ [V]")
ax.set_ylabel(r"$V^\textrm{in}_2$ [V]")

plt.colorbar(im)

plt.savefig(
    join(pathToSimFolder, "fitness_2D_normed.png"), bbox_inches="tight", dpi=300
)
# plt.show()
plt.close(fig)

if voltageScanPoints != 2:
    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.errorbar(
        np.arange(4),
        currents2x2.flatten(),
        yerr=currentUncert2x2.flatten(),
        marker="x",
        color="k",
        linestyle="",
    )

    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels([r"[0,0]", r"[1,0]", r"[0,1]", r"[1,1]"])

    plt.savefig(join(pathToSimFolder, "fitness_1D.png"), bbox_inches="tight", dpi=300)

    ax.set_xlabel("gate input")
    ax.set_ylabel(r"$I^{\textrm{out}} [e \nu_0]$")
    # plt.show()
    plt.close(fig)
