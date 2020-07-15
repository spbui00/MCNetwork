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
    print(inp[fileNumber - 1])
    # for fileNumber in [1]:

    data = np.genfromtxt(
        join(pathToSimFolder, f"swapTrackFile{fileNumber}.txt"),
        delimiter=";",
        dtype=int,
    )

    trajectoriesSortedByStartEnd = [
        [[] for j in range(len(electrodes))] for i in range(len(electrodes))
    ]
    trajectories = []
    hops = 20000
    IDs = {}
    hitID = 0
    for i in range(hops):
        hoppingSite1 = data[i, 0]
        hoppingSite2 = data[i, 1]
        # print("hoppingSite1",hoppingSite1,"hoppingSite2",hoppingSite2)
        if hoppingSite1 in IDs:
            ID = IDs[hoppingSite1]
            del IDs[hoppingSite1]
            # print("found ID",ID)
        else:
            ID = hitID
            hitID += 1
            trajectories.append([])
            # print("new ID", ID)

        if hoppingSite2 < parameters["acceptorNumber"]:
            IDs[hoppingSite2] = ID

        trajectories[ID].append([hoppingSite1, hoppingSite2])

    # sort trajectories
    for i in range(len(trajectories)):
        if trajectories[i][0][0] >= parameters["acceptorNumber"]:
            if trajectories[i][-1][1] >= parameters["acceptorNumber"]:
                trajectoriesSortedByStartEnd[
                    trajectories[i][0][0] - int(parameters["acceptorNumber"])
                ][trajectories[i][-1][1] - int(parameters["acceptorNumber"])].append(
                    trajectories[i]
                )
                # print(trajectories[i][0][0], trajectories[i][-1][1])

    for k in range(len(electrodes)):

        fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

        electodePlotWidth = 8
        for i in range(len(electrodes)):
            if i == parameters["outputElectrode"]:
                col = "blue"
            elif i == parameters["inputElectrode1"]:
                if fileNumber in [3, 4]:
                    col = "red"
                else:
                    col = "rosybrown"
            elif i == parameters["inputElectrode2"]:
                if fileNumber in [2, 4]:
                    col = "red"
                else:
                    col = "rosybrown"
            else:
                col = "green"

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
                        fc=col,
                        ec=col,
                        zorder=-1,
                    )
                )
            elif parameters["geometry"] == "circle":
                electrodeWidth = (
                    parameters["electrodeWidth"]
                    / (parameters["radius"] * 2 * np.pi)
                    * 360
                )  # in degrees
                ax.add_artist(
                    Wedge(
                        (0, 0),
                        parameters["radius"] + electodePlotWidth / 2,
                        electrodes[i][0] - electrodeWidth / 2,
                        electrodes[i][0] + electrodeWidth / 2,
                        width=electodePlotWidth,
                        fc=col,
                        ec=col,
                        zorder=-1,
                    )
                )

        ax.scatter(acceptorPos[:, 0], acceptorPos[:, 1], c="k", marker=".", s=20)
        ax.scatter(donorPos[:, 0], donorPos[:, 1], c="k", marker="x", s=20)

        for l in range(len(electrodes)):
            trajectories = trajectoriesSortedByStartEnd[k][l]

            for i in range(len(trajectories)):
                for j in range(len(trajectories[i])):

                    hoppingSite1 = trajectories[i][j][0]
                    hoppingSite2 = trajectories[i][j][1]

                    if hoppingSite1 >= parameters["acceptorNumber"]:
                        x1, y1 = (
                            electrodePositions[
                                hoppingSite1 - int(parameters["acceptorNumber"])
                            ][0],
                            electrodePositions[
                                hoppingSite1 - int(parameters["acceptorNumber"])
                            ][1],
                        )
                    else:
                        x1, y1 = (
                            acceptorPos[hoppingSite1, 0],
                            acceptorPos[hoppingSite1, 1],
                        )

                    if hoppingSite2 >= parameters["acceptorNumber"]:
                        x2, y2 = (
                            electrodePositions[
                                hoppingSite2 - int(parameters["acceptorNumber"])
                            ][0],
                            electrodePositions[
                                hoppingSite2 - int(parameters["acceptorNumber"])
                            ][1],
                        )
                    else:
                        x2, y2 = (
                            acceptorPos[hoppingSite2, 0],
                            acceptorPos[hoppingSite2, 1],
                        )

                    # ax.plot([x1,x2],[y1,y2],"-",alpha=0.05,color="k",linewidth=2)
                    ax.plot(
                        [x1, x2],
                        [y1, y2],
                        "-",
                        alpha=0.05,
                        color=color(l, len(electrodes)),
                        linewidth=2,
                    )

                    # if currentRatio>0.5:
                    #     ax.arrow((x2+x1)/2,(y2+y1)/2,(x2-x1)*0.001,(y2-y1)*0.001,color=colorMaker(abs(currentRatio-0.5)*2),ec=None,alpha=absBins[i,j],linewidth=0,head_width=(currentRatio-0.5)*20)

        ax.axis("off")
        if parameters["geometry"] == "circle":
            ax.add_artist(
                plt.Circle((0, 0), parameters["radius"], fc="none", ec="k", zorder=-2)
            )
        elif parameters["geometry"] == "rect":
            ax.add_artist(
                plt.Rectangle(
                    (0, 0),
                    parameters["lenX"],
                    parameters["lenY"],
                    fc="none",
                    ec="k",
                    zorder=-2,
                )
            )

        if parameters["geometry"] == "rect":
            ax.set_xlim(
                -electodePlotWidth / 2, parameters["lenX"] + electodePlotWidth / 2
            )
            ax.set_ylim(
                -electodePlotWidth / 2, parameters["lenY"] + electodePlotWidth / 2
            )
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

        plt.savefig(
            join(pathToSimFolder, f"trajectory_fromEl_{k}_{inp[fileNumber-1]}.png"),
            bbox_inches="tight",
            dpi=300,
        )

        # plt.show()
        plt.close(fig)

    for k in range(len(electrodes)):

        fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

        electodePlotWidth = 8
        for i in range(len(electrodes)):
            if i == parameters["outputElectrode"]:
                col = "blue"
            elif i == parameters["inputElectrode1"]:
                if fileNumber in [3, 4]:
                    col = "red"
                else:
                    col = "rosybrown"
            elif i == parameters["inputElectrode2"]:
                if fileNumber in [2, 4]:
                    col = "red"
                else:
                    col = "rosybrown"
            else:
                col = "green"

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
                        fc=col,
                        ec=col,
                        zorder=-1,
                    )
                )
            elif parameters["geometry"] == "circle":
                electrodeWidth = (
                    parameters["electrodeWidth"]
                    / (parameters["radius"] * 2 * np.pi)
                    * 360
                )  # in degrees
                ax.add_artist(
                    Wedge(
                        (0, 0),
                        parameters["radius"] + electodePlotWidth / 2,
                        electrodes[i][0] - electrodeWidth / 2,
                        electrodes[i][0] + electrodeWidth / 2,
                        width=electodePlotWidth,
                        fc=col,
                        ec=col,
                        zorder=-1,
                    )
                )

        ax.scatter(acceptorPos[:, 0], acceptorPos[:, 1], c="k", marker=".", s=20)
        ax.scatter(donorPos[:, 0], donorPos[:, 1], c="k", marker="x", s=20)

        for l in range(len(electrodes)):
            trajectories = trajectoriesSortedByStartEnd[l][k]

            for i in range(len(trajectories)):
                for j in range(len(trajectories[i])):

                    hoppingSite1 = trajectories[i][j][0]
                    hoppingSite2 = trajectories[i][j][1]

                    if hoppingSite1 >= parameters["acceptorNumber"]:
                        x1, y1 = (
                            electrodePositions[
                                hoppingSite1 - int(parameters["acceptorNumber"])
                            ][0],
                            electrodePositions[
                                hoppingSite1 - int(parameters["acceptorNumber"])
                            ][1],
                        )
                    else:
                        x1, y1 = (
                            acceptorPos[hoppingSite1, 0],
                            acceptorPos[hoppingSite1, 1],
                        )

                    if hoppingSite2 >= parameters["acceptorNumber"]:
                        x2, y2 = (
                            electrodePositions[
                                hoppingSite2 - int(parameters["acceptorNumber"])
                            ][0],
                            electrodePositions[
                                hoppingSite2 - int(parameters["acceptorNumber"])
                            ][1],
                        )
                    else:
                        x2, y2 = (
                            acceptorPos[hoppingSite2, 0],
                            acceptorPos[hoppingSite2, 1],
                        )

                    # ax.plot([x1,x2],[y1,y2],"-",alpha=0.05,color="k",linewidth=2)
                    ax.plot(
                        [x1, x2],
                        [y1, y2],
                        "-",
                        alpha=0.05,
                        color=color(l, len(electrodes)),
                        linewidth=2,
                    )

                    # if currentRatio>0.5:
                    #     ax.arrow((x2+x1)/2,(y2+y1)/2,(x2-x1)*0.001,(y2-y1)*0.001,color=colorMaker(abs(currentRatio-0.5)*2),ec=None,alpha=absBins[i,j],linewidth=0,head_width=(currentRatio-0.5)*20)

        ax.axis("off")
        if parameters["geometry"] == "circle":
            ax.add_artist(
                plt.Circle((0, 0), parameters["radius"], fc="none", ec="k", zorder=-2)
            )
        elif parameters["geometry"] == "rect":
            ax.add_artist(
                plt.Rectangle(
                    (0, 0),
                    parameters["lenX"],
                    parameters["lenY"],
                    fc="none",
                    ec="k",
                    zorder=-2,
                )
            )

        if parameters["geometry"] == "rect":
            ax.set_xlim(
                -electodePlotWidth / 2, parameters["lenX"] + electodePlotWidth / 2
            )
            ax.set_ylim(
                -electodePlotWidth / 2, parameters["lenY"] + electodePlotWidth / 2
            )
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

        plt.savefig(
            join(pathToSimFolder, f"trajectory_toEl_{k}_{inp[fileNumber-1]}.png"),
            bbox_inches="tight",
            dpi=300,
        )

        # plt.show()
        plt.close(fig)
