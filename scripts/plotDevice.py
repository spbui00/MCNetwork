#!/usr/bin/python3

from sys import argv
from tools import *
from os.path import join

import numpy as np
from matplotlib.patches import Wedge

if len(argv) > 1:
    pathToSimFolder = argv[1]
else:
    pathToSimFolder = "../data/"

parameters, electrodes = readParameters(pathToSimFolder)


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
print("plotting device")


import matplotlib.pylab as plt


fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))


electodePlotWidth = 8
for i in range(len(electrodes)):
    if i == parameters["outputElectrode"]:
        color = "blue"
    elif i == parameters["inputElectrode1"]:
        color = "red"
    elif i == parameters["inputElectrode2"]:
        color = "red"
    else:
        color = "green"

    if parameters["geometry"] == "rect":
        if electrodes[i][1] == 0:
            angle = 0
            xy = (
                0 - electodePlotWidth / 2,
                electrodes[i][0] * parameters["lenY"]
                - parameters["electrodeWidth"] / 2,
            )
        elif electrodes[i][1] == 1:
            angle = 0
            xy = (
                parameters["lenX"] - electodePlotWidth / 2,
                electrodes[i][0] * parameters["lenY"]
                - parameters["electrodeWidth"] / 2,
            )
        elif electrodes[i][1] == 2:
            angle = 90
            xy = (
                electrodes[i][0] * parameters["lenX"]
                + parameters["electrodeWidth"] / 2,
                0 - electodePlotWidth / 2,
            )
        elif electrodes[i][1] == 3:
            angle = 90
            xy = (
                electrodes[i][0] * parameters["lenX"]
                + parameters["electrodeWidth"] / 2,
                parameters["lenY"] - electodePlotWidth / 2,
            )
        ax.add_artist(
            plt.Rectangle(
                xy,
                electodePlotWidth,
                parameters["electrodeWidth"],
                angle=angle,
                fc=color,
                ec=color,
                zorder=-1,
            )
        )
    elif parameters["geometry"] == "circle":
        electrodeWidth = (
            parameters["electrodeWidth"] / (parameters["radius"] * 2 * np.pi) * 360
        )  # in degrees
        ax.add_artist(
            Wedge(
                (0, 0),
                parameters["radius"] + electodePlotWidth / 2,
                electrodes[i][0] - electrodeWidth / 2,
                electrodes[i][0] + electrodeWidth / 2,
                width=electodePlotWidth,
                fc=color,
                ec=color,
                zorder=-1,
            )
        )


ax.scatter(acceptorPos[:, 0], acceptorPos[:, 1], c="k", marker=".", s=20)
ax.scatter(donorPos[:, 0], donorPos[:, 1], c="k", marker="x", s=20)

ax.axis("off")

if parameters["geometry"] == "circle":
    ax.add_artist(
        plt.Circle((0, 0), parameters["radius"], fc="none", ec="k", zorder=-2)
    )
elif parameters["geometry"] == "rect":
    ax.add_artist(
        plt.Rectangle(
            (0, 0), parameters["lenX"], parameters["lenY"], fc="none", ec="k", zorder=-2
        )
    )


if parameters["geometry"] == "rect":
    ax.set_xlim(-electodePlotWidth / 2, parameters["lenX"] + electodePlotWidth / 2)
    ax.set_ylim(-electodePlotWidth / 2, parameters["lenY"] + electodePlotWidth / 2)
elif parameters["geometry"] == "circle":
    ax.set_xlim(
        -parameters["radius"] - electodePlotWidth,
        parameters["radius"] + electodePlotWidth,
    )
    ax.set_ylim(
        -parameters["radius"] - electodePlotWidth,
        parameters["radius"] + electodePlotWidth,
    )

ax.set_aspect("equal")


plt.savefig(join(pathToSimFolder, "plotDevice.png"), bbox_inches="tight", dpi=300)
# plt.show()
plt.close()
fig = None
# print(parameters)
