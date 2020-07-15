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
            voltages = np.array(dataFile["/voltages"][:])
            optEnergy = np.array(dataFile["/optEnergy"][:])
            while True:
                try:
                    generations = np.array(dataFile["/generation"][:])
                    mode = "genetic"
                    break
                except KeyError:
                    pass
                try:
                    basinAccepted = np.array(dataFile["/basinAccepted"][:], dtype=int)
                    accepted = basinAccepted.astype(bool)
                    notAccepted = np.invert(accepted)
                    mode = "basinHop"
                    break
                except KeyError:
                    pass
                mode = "MC"
                try:
                    accepted = np.array(dataFile["/accepted"][:], dtype=bool)
                    notAccepted = np.invert(accepted)
                except KeyError:
                    accepted = np.ones(
                        optEnergy.shape, dtype=bool
                    )  # support for deprecated version
                    notAccepted = np.invert(accepted)
                break
        break
    except OSError as e:
        if "No such file" in repr(e):
            raise e
        else:
            print(f"could not open file. try number {fileOpenTries}")
            sleep(1)


cotrolElectrodeIndices = list(range(0, len(electrodes)))
cotrolElectrodeIndices.remove(parameters["outputElectrode"])
cotrolElectrodeIndices.remove(parameters["inputElectrode1"])
cotrolElectrodeIndices.remove(parameters["inputElectrode2"])

controlVoltages = voltages[:, cotrolElectrodeIndices]


if mode == "MC":

    distance = 0
    meanRange = 1000
    displace = []

    for i in range(int(distance + meanRange / 2), controlVoltages.shape[0]):
        mean = np.mean(
            controlVoltages[
                int(i - distance - meanRange / 2) : int(i - distance + meanRange / 2), :
            ],
            axis=0,
        )
        # displace.append(np.sqrt(np.sum((controlVoltages[i])**2)))
        displace.append(np.sqrt(np.sum((mean - controlVoltages[i]) ** 2)))

    MSD = np.sum((controlVoltages[0] - controlVoltages[:]) ** 2, axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(range(len(MSD)), MSD, "r-", label="MSD")
    ax2 = ax.twinx()
    ax2.plot(
        range(int(distance + meanRange / 2), controlVoltages.shape[0]),
        displace,
        "k-",
        label="displacement",
    )

    ax.legend()
    ax2.legend()

    ax.set_xlabel("step")
    ax.set_ylabel("displacement")

    plt.savefig(join(pathToSimFolder, "displacement.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(
        np.arange(optEnergy.shape[0])[notAccepted[:, 0]],
        optEnergy[notAccepted],
        ".",
        ms=1,
        color="darkred",
        label="not accepted",
        zorder=10,
    )
    ax.plot(
        np.arange(optEnergy.shape[0])[accepted[:, 0]],
        optEnergy[accepted],
        ".",
        ms=1,
        color="darkgreen",
        label="accepted",
        zorder=10,
    )

    # ax.set_xlim(-0.15,0.65)
    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")

    ax.legend()

    plt.savefig(join(pathToSimFolder, "convergence.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(
        np.arange(optEnergy.shape[0])[notAccepted[:, 0]],
        optEnergy[notAccepted],
        ".",
        ms=1,
        color="darkred",
        label="not accepted",
        zorder=10,
    )
    ax.plot(
        np.arange(optEnergy.shape[0])[accepted[:, 0]],
        optEnergy[accepted],
        ".",
        ms=1,
        color="darkgreen",
        label="accepted",
        zorder=10,
    )

    ax2 = ax.twinx()
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    ax2.plot(
        range(int(distance + meanRange / 2), controlVoltages.shape[0]),
        displace,
        "k-",
        label="displacement",
    )

    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")
    ax2.set_ylabel("displacement")

    # ax.legend([line],[line.get_label()])
    # ax2.legend()

    plt.savefig(
        join(pathToSimFolder, "convergence_displacement.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)


###############################


if mode == "genetic":

    distance = 0
    meanRange = 1000
    displace = []

    for i in range(int(distance + meanRange / 2), controlVoltages.shape[0]):
        mean = np.mean(
            controlVoltages[
                int(i - distance - meanRange / 2) : int(i - distance + meanRange / 2), :
            ],
            axis=0,
        )
        # displace.append(np.sqrt(np.sum((controlVoltages[i])**2)))
        displace.append(np.sqrt(np.sum((mean - controlVoltages[i]) ** 2)))

    MSD = np.sum((controlVoltages[0] - controlVoltages[:]) ** 2, axis=1)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(range(len(MSD)), MSD, "r-", label="MSD")
    ax2 = ax.twinx()
    ax2.plot(
        range(int(distance + meanRange / 2), controlVoltages.shape[0]),
        displace,
        "k-",
        label="displacement",
    )

    ax.legend()
    ax2.legend()

    ax.set_xlabel("step")
    ax.set_ylabel("displacement")

    plt.savefig(join(pathToSimFolder, "displacement.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)

    genBest = np.empty(optEnergy.shape)
    for i in range(int(optEnergy.shape[0] / 25)):
        genBest[i * 25 : (i + 1) * 25] = max(optEnergy[i * 25 : (i + 1) * 25])

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(optEnergy, ".", ms=1, color="darkgreen", label="all")
    ax.plot(genBest, color="darkblue", label="gen best")

    # ax.set_xlim(-0.15,0.65)
    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")

    ax.legend()

    plt.savefig(join(pathToSimFolder, "convergence.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(optEnergy, ".", ms=1, color="darkgreen", label="all")
    ax.plot(genBest, color="darkblue", label="gen best")

    ax2 = ax.twinx()
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    ax2.plot(
        range(int(distance + meanRange / 2), controlVoltages.shape[0]),
        displace,
        "k-",
        label="displacement",
    )

    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")
    ax2.set_ylabel("displacement")

    # ax.legend([line],[line.get_label()])
    # ax2.legend()

    plt.savefig(
        join(pathToSimFolder, "convergence_displacement.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)


###############################
if mode == "basinHop":

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(
        np.arange(optEnergy.shape[0])[notAccepted[:, 0]],
        optEnergy[notAccepted],
        ".",
        ms=1,
        color="darkred",
        label="not accepted",
        zorder=10,
    )
    ax.plot(
        np.arange(optEnergy.shape[0])[accepted[:, 0]],
        optEnergy[accepted],
        ".",
        ms=1,
        color="darkgreen",
        label="accepted",
        zorder=10,
    )

    buff = np.where(basinAccepted[:, 0] == 2)[0]
    basinChanges = np.array([buff, np.zeros(buff.shape)], dtype=int)
    buff = np.where(basinAccepted[:, 0] == 3)[0]
    basinChanges = np.append(
        basinChanges, np.array([buff, np.ones(buff.shape)], dtype=int), axis=1
    )
    basinChanges = basinChanges[:, np.argsort(basinChanges[0])]

    if basinChanges.shape[1] > 0:
        for i in range(basinChanges.shape[1]):
            if basinChanges[1, i]:
                ax.axvline(basinChanges[0, i], color="darkgreen", zorder=-1)
            else:
                ax.axvline(basinChanges[0, i], color="darkred", zorder=-1)

        ax.plot(
            np.arange(0, basinChanges[0, 0]),
            np.maximum.accumulate(optEnergy[: basinChanges[0, 0]]),
            color="darkblue",
            label="basin best",
        )
        for i in range(1, basinChanges.shape[1]):
            ax.plot(
                np.arange(basinChanges[0, i - 1], basinChanges[0, i]),
                np.maximum.accumulate(
                    optEnergy[basinChanges[0, i - 1] : basinChanges[0, i]]
                ),
                color="darkblue",
            )
        ax.plot(
            np.arange(basinChanges[0, -1], len(optEnergy)),
            np.maximum.accumulate(optEnergy[basinChanges[0, -1] :]),
            color="darkblue",
        )

    # ax.set_xlim(-0.15,0.65)
    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")

    ax.legend()

    plt.savefig(join(pathToSimFolder, "convergence.png"), bbox_inches="tight", dpi=300)
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(np.maximum.accumulate(optEnergy), color="darkorange", label="best")

    ax.plot(
        np.arange(optEnergy.shape[0])[notAccepted[:, 0]],
        optEnergy[notAccepted],
        ".",
        ms=1,
        color="darkred",
        label="not accepted",
        zorder=10,
    )
    ax.plot(
        np.arange(optEnergy.shape[0])[accepted[:, 0]],
        optEnergy[accepted],
        ".",
        ms=1,
        color="darkgreen",
        label="accepted",
        zorder=10,
    )

    if basinChanges.shape[1] > 0:
        for i in range(basinChanges.shape[1]):
            if basinChanges[1, i]:
                ax.axvline(basinChanges[0, i], color="darkgreen", zorder=-1)
            else:
                ax.axvline(basinChanges[0, i], color="darkred", zorder=-1)

        ax.plot(
            np.arange(0, basinChanges[0, 0]),
            np.maximum.accumulate(optEnergy[: basinChanges[0, 0]]),
            color="darkblue",
            label="basin best",
        )
        for i in range(1, basinChanges.shape[1]):
            ax.plot(
                np.arange(basinChanges[0, i - 1], basinChanges[0, i]),
                np.maximum.accumulate(
                    optEnergy[basinChanges[0, i - 1] : basinChanges[0, i]]
                ),
                color="darkblue",
            )
        ax.plot(
            np.arange(basinChanges[0, -1], len(optEnergy)),
            np.maximum.accumulate(optEnergy[basinChanges[0, -1] :]),
            color="darkblue",
        )

    ax2 = ax.twinx()
    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    # calc last basin best
    basinBestIdx = np.argmax(optEnergy[0 : basinChanges[0, 0]])
    basinBestVoltages = controlVoltages[basinBestIdx]

    # ax2.plot(np.arange(0,basinChanges[0,0]), np.sqrt(np.sum((controlVoltages[0:basinChanges[0,0]] - basinBestVoltages)**2, axis = 1  )) ,color="darkblue")
    for i in range(1, basinChanges.shape[1]):
        ax2.plot(
            np.arange(basinChanges[0, i - 1], basinChanges[0, i]),
            np.sqrt(
                np.sum(
                    (
                        controlVoltages[basinChanges[0, i - 1] : basinChanges[0, i]]
                        - basinBestVoltages
                    )
                    ** 2,
                    axis=1,
                )
            ),
            color="k",
        )
        # calc last basin best
        if basinChanges[1, i]:
            basinBestIdx = (
                np.argmax(optEnergy[basinChanges[0, i - 1] : basinChanges[0, i]])
                + basinChanges[0, i - 1]
            )
            basinBestVoltages = controlVoltages[basinBestIdx]
    ax2.plot(
        np.arange(basinChanges[0, -1], len(optEnergy)),
        np.sqrt(
            np.sum(
                (controlVoltages[basinChanges[0, -1] :] - basinBestVoltages) ** 2,
                axis=1,
            )
        ),
        color="k",
    )

    ax.set_ylim(0.15, 1.05)

    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$\mathcal{F}$")
    ax2.set_ylabel("dist")

    # ax.legend([line],[line.get_label()])
    # ax2.legend()

    plt.savefig(
        join(pathToSimFolder, "convergence_dist.png"), bbox_inches="tight", dpi=300
    )
    # plt.show()
    plt.close(fig)
