#!/usr/bin/python3
"""
copys best optEnergy voltages to input file
"""
from tools import *
from sys import argv
from os.path import join

import h5py
import numpy as np

import shutil

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
        break
    except OSError as e:
        if "No such file" in repr(e):
            raise e
        else:
            print(f"could not open file. try number {fileOpenTries}")
            sleep(1)


best = np.argmax(optEnergy, axis=0)
# best=np.argmax(fitness,axis=0)
# best=[0]
print("best at", best, "optEnergy: ", optEnergy[best])
print(voltages[best])


i = 0
with open(join(pathToSimFolder, "in.txt")) as oldFile:
    with open(join(pathToSimFolder, "in2.txt"), "w") as newFile:
        for rawLine in oldFile:
            # print(rawLine)
            line = sub(r" #.*", "", rawLine)
            line = sub(r"#.*", "", line)

            splitted = line.split(" ")
            # print(splitted)
            if splitted[0] == "electrode":
                try:
                    geom = parameters["geometry"]
                except KeyError:
                    geom = "rect"
                if geom == "rect":
                    newFile.write(
                        str(splitted[0])
                        + " "
                        + str(splitted[1])
                        + " "
                        + str(splitted[2])
                        + " "
                        + str(voltages[best][0, i])
                        + "\n"
                    )
                elif geom == "circle":
                    newFile.write(
                        str(splitted[0])
                        + " "
                        + str(splitted[1])
                        + " "
                        + str(voltages[best][0, i])
                        + "\n"
                    )
                i += 1
            else:
                newFile.write(rawLine)

shutil.move(join(pathToSimFolder, "in2.txt"), join(pathToSimFolder, "in.txt"))
