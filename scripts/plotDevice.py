#!/usr/bin/python3

from sys import argv
from tools import readParameters
from os.path import join

import numpy as np

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)


acceptorPos=np.zeros((int(parameters["acceptorNumber"]),2))
donorPos   =np.zeros((int(parameters["acceptorNumber"]*parameters["compensationFactor"]),2))
with open(join(pathToSimFolder,"device.txt")) as deviceFile:
    line=next(deviceFile)
    for i in range(acceptorPos.shape[0]):
        acceptorPos[i]=next(deviceFile).split(" ")
    line=next(deviceFile)
    line=next(deviceFile)
    for i in range(donorPos.shape[0]):
        donorPos[i]=next(deviceFile).split(" ")

print(acceptorPos)
print(donorPos)


import matplotlib.pylab as plt


fig, ax = plt.subplots(1,1)

for i in range(len(electrodes)):
    if   i == parameters["outputElectrode"]: color= "blue"
    elif i == parameters["inputElectrode1"]: color= "red"
    elif i == parameters["inputElectrode2"]: color= "red"
    else:                                    color= "green"
    
    if electrodes[i][1]==0: ax.scatter(0                              ,electrodes[i][0]*parameters["lenY"],c=color,marker=".",s=400)
    if electrodes[i][1]==1: ax.scatter(parameters["lenX"]             ,electrodes[i][0]*parameters["lenY"],c=color,marker=".",s=400)
    if electrodes[i][1]==2: ax.scatter(electrodes[i][0]*parameters["lenX"],0                              ,c=color,marker=".",s=400)
    if electrodes[i][1]==3: ax.scatter(electrodes[i][0]*parameters["lenX"],parameters["lenY"]             ,c=color,marker=".",s=400)


for i in range(acceptorPos.shape[0]):
    ax.text(acceptorPos[i,0],acceptorPos[i,1], f"{i}",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=10, color='red',
            transform=ax.transData)
    ax.add_artist(plt.Circle((acceptorPos[i,0],acceptorPos[i,1]),8,fc='none',ec="r"))

for i in range(donorPos.shape[0]):
    ax.text(donorPos[i,0],donorPos[i,1], f"{i}",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=10, color='blue',
            transform=ax.transData)
    ax.add_artist(plt.Circle((donorPos[i,0],donorPos[i,1]),8,fc='none',ec="b"))


ax.set_xlim(0,parameters["lenX"])
ax.set_ylim(0,parameters["lenY"])

ax.set_aspect('equal')

ax.set_xticks([], [])
ax.set_yticks([], [])

plt.show()

# print(parameters)
