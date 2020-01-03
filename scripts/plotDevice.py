#!/usr/bin/python3

from sys import argv
from tools import readParameters

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)




import matplotlib.pylab as plt


fig, ax = plt.subplots(1,1)

for electrode in electrodes:
    ax.scatter()

ax.set


print(parameters)
