#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
from matplotlib.patches import Wedge

import numpy as np

if len(argv) > 1:
    pathToSimFolder = argv[1]
else:
    pathToSimFolder = "../data/"

parameters, electrodes = readParameters(pathToSimFolder)
electrodeNumber = len(electrodes)

acceptorPos = np.zeros((int(parameters["acceptorNumber"]), 2))
try:
    donorPos = np.zeros((int(parameters["donorNumber"]), 2))
except KeyError:
    donorPos = np.zeros(
        (int(parameters["acceptorNumber"] * parameters["compensationFactor"]), 2)
    )

with open(join(pathToSimFolder, "device.txt")) as deviceFile:
    line = next(deviceFile)
    line = next(deviceFile)
    for i in range(acceptorPos.shape[0]):
        acceptorPos[i] = next(deviceFile).split(" ")
    line = next(deviceFile)
    line = next(deviceFile)
    for i in range(donorPos.shape[0]):
        donorPos[i] = next(deviceFile).split(" ")

# print(acceptorPos)
# print(donorPos)

electrodePositions = np.empty((len(electrodes), 2))
for i in range(len(electrodes)):
    if parameters["geometry"] == "rect":
        if electrodes[i][1] == 0:
            electrodePositions[i] = [0, electrodes[i][0] * parameters["lenY"]]
        if electrodes[i][1] == 1:
            electrodePositions[i] = [
                parameters["lenX"],
                electrodes[i][0] * parameters["lenY"],
            ]
        if electrodes[i][1] == 2:
            electrodePositions[i] = [electrodes[i][0] * parameters["lenX"], 0]
        if electrodes[i][1] == 3:
            electrodePositions[i] = [
                electrodes[i][0] * parameters["lenX"],
                parameters["lenY"],
            ]
    elif parameters["geometry"] == "circle":
        electrodePositions[i] = [
            parameters["radius"] * np.cos(electrodes[i][0] / 360 * 2 * np.pi),
            parameters["radius"] * np.sin(electrodes[i][0] / 360 * 2 * np.pi),
        ]

# print(electrodePositions)


def colorMaker(x):
    from matplotlib import colors
    from scipy.interpolate import interp1d

    cols = ["darkred", "darkgreen"]

    rgbaData = np.array([colors.to_rgba(c) for c in cols])
    rInterpolater = interp1d(np.linspace(0, 1, len(cols)), rgbaData[:, 0])
    gInterpolater = interp1d(np.linspace(0, 1, len(cols)), rgbaData[:, 1])
    bInterpolater = interp1d(np.linspace(0, 1, len(cols)), rgbaData[:, 2])
    return np.array([rInterpolater(x), gInterpolater(x), bInterpolater(x), 1])


inp = ["0_0", "0_1", "1_0", "1_1"]
for fileNumber in [1, 2, 3, 4]:
    # for fileNumber in [1]:

    data = np.genfromtxt(
        join(pathToSimFolder, f"swapTrackFile{fileNumber}.txt"),
        delimiter=";",
        dtype=int,
    )

    maxIndex = np.max(data)
    added = (maxIndex + 1) * data[:, 0] + data[:, 1]

    bins = np.bincount(added)
    bins.resize(maxIndex + 1, maxIndex + 1)

    absBins = bins + bins.T
    absBins = absBins / np.max(absBins)

    distances = []
    swapps = []

    for i in range(bins.shape[0]):
        if i >= parameters["acceptorNumber"]:
            x1, y1 = (
                electrodePositions[i - int(parameters["acceptorNumber"])][0],
                electrodePositions[i - int(parameters["acceptorNumber"])][1],
            )
        else:
            x1, y1 = acceptorPos[i, 0], acceptorPos[i, 1]

        for j in range(i):
            if j >= parameters["acceptorNumber"]:
                x2, y2 = (
                    electrodePositions[j - int(parameters["acceptorNumber"])][0],
                    electrodePositions[j - int(parameters["acceptorNumber"])][1],
                )
            else:
                x2, y2 = acceptorPos[j, 0], acceptorPos[j, 1]

            # ax.plot([x1,x2],[y1,y2],"k-",alpha=bins[i,j])
            if (bins[i, j] + bins[j, i]) != 0:
                currentRatio = bins[i, j] / (bins[i, j] + bins[j, i])

                distances.append(np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2))
                swapps.append(absBins[i, j])

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        distances,
        weights=swapps,
        bins=50,
        density=True,
        histtype="step",
        color="k",
        range=[0, 100],
    )

    ax.set_xlabel(r"d [nm]")
    ax.set_ylabel(r"P(swapp)")

    plt.savefig(
        join(pathToSimFolder, f"swapHist_{inp[fileNumber-1]}.png"),
        bbox_inches="tight",
        dpi=300,
    )

    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        distances,
        weights=swapps,
        bins=50,
        density=True,
        histtype="step",
        color="k",
        range=[0, 100],
        log=True,
    )

    ax.set_xlabel(r"d [nm]")
    ax.set_ylabel(r"P(swapp)")

    plt.savefig(
        join(pathToSimFolder, f"swapHist_{inp[fileNumber-1]}_log.png"),
        bbox_inches="tight",
        dpi=300,
    )

    # plt.show()
    plt.close(fig)
