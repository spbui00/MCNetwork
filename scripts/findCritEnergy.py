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


for distance in [25, 50, 100, 200, 500, 1000]:
    runs = 40

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    xRange = np.arange(0.4, 1, 0.0001)
    data = np.ma.array(np.zeros((runs, *xRange.shape)), mask=True)

    for r in range(runs):
        print(r + 1)
        fileOpenTries = 0
        while fileOpenTries < 50:
            fileOpenTries += 1
            try:
                with h5py.File(
                    join(pathToSimFolder, f"run_{r+1}/data.hdf5"), "r"
                ) as dataFile:
                    optEnergy = np.array(dataFile["/optEnergy"][20:, 0])
                break
            except OSError as e:
                if "No such file" in repr(e):
                    raise e
                else:
                    print(f"could not open file. try number {fileOpenTries}")
                    sleep(1)

        accumulated = np.maximum.accumulate(optEnergy)

        dataX = []
        dataY = []

        currentValue = 0
        shiftedValue = 0

        N = accumulated.shape[0]
        for i in range(N):
            if (
                accumulated[i] != currentValue
                or accumulated[min(i + distance, N - 1)] != shiftedValue
            ):
                currentValue = accumulated[i]
                shiftedValue = accumulated[min(i + distance, N - 1)]

                dataX.append(currentValue)
                dataY.append(shiftedValue)

        """
        fig2, ax2=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


        ax2.plot(dataX, dataY,"k-")


        ax2.set_xlim(0.45,1)
        ax2.set_ylim(0.45,1)

        ax2.set_xlabel(r"$\mathcal{F}_{i}$")
        ax2.set_ylabel(rf"$\mathcal{{F}}_{{i+{distance} }}$")


        plt.savefig(join(pathToSimFolder,f"run_{r+1}/critFitness_{r+1}.png"),bbox_inches="tight",dpi=300)
        # plt.show()
        plt.close(fig)
        """

        dataX = np.array(dataX)
        dataY = np.array(dataY)

        ax.plot(dataX, dataY, "-", color=color(r, runs))

        # stuff for mean plot
        indices = np.digitize(xRange, dataX)
        firstOut = np.where(indices == dataX.shape[0])[0][0]
        indices = indices[:firstOut]
        data[r, :firstOut] = dataY[indices]

    ax.plot([0.4, 1], [0.4, 1], "k-", lw=1)

    ax.set_xlim(0.4, 1)
    ax.set_ylim(0.4, 1)

    ax.set_xlabel(r"$\mathcal{F}_{i}$")
    ax.set_ylabel(rf"$\mathcal{{F}}_{{i+{distance} }}$")

    plt.savefig(
        join(pathToSimFolder, f"critFitness_{distance}.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    mean = np.mean(data, axis=0)
    uncert = StdMean(data, axis=0)
    ax.plot(xRange, mean, "k-")
    ax.fill_between(xRange, mean - uncert, mean + uncert, facecolor="grey")

    ax.plot([0.4, 1], [0.4, 1], "k-", lw=1)

    ax.set_xlim(0.4, 1)
    ax.set_ylim(0.4, 1)

    ax.set_xlabel(r"$\mathcal{F}_{i}$")
    ax.set_ylabel(rf"$\mathcal{{F}}_{{i+{distance} }}$")

    plt.savefig(
        join(pathToSimFolder, f"critFitness_mean_{distance}.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)
