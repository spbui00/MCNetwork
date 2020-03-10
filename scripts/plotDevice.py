#!/usr/bin/python3

from sys import argv
from tools import *
from os.path import join

import numpy as np
from matplotlib.patches import Wedge, Rectangle

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)


acceptorPos = np.zeros((int(parameters["acceptorNumber"]),2))
try:
    donorPos    = np.zeros((int(parameters["donorNumber"]),2))
except KeyError:
    donorPos    = np.zeros((int(parameters["acceptorNumber"]*parameters["compensationFactor"]),2))
    
with open(join(pathToSimFolder,"device.txt")) as deviceFile:
    line=next(deviceFile)
    line=next(deviceFile)
    for i in range(acceptorPos.shape[0]):
        acceptorPos[i]=next(deviceFile).split(" ")
    line=next(deviceFile)
    line=next(deviceFile)
    for i in range(donorPos.shape[0]):
        donorPos[i]=next(deviceFile).split(" ")

print(acceptorPos)
print(donorPos)


electrodePositions=np.empty((len(electrodes),2))
for i in range(len(electrodes)):
    if parameters["geometry"] == "rect":
        if electrodes[i][1]==0: electrodePositions[i] = [0                              ,electrodes[i][0]*parameters["lenY"]]
        if electrodes[i][1]==1: electrodePositions[i] = [parameters["lenX"]             ,electrodes[i][0]*parameters["lenY"]]
        if electrodes[i][1]==2: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],0                              ]
        if electrodes[i][1]==3: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],parameters["lenY"]             ]
    elif parameters["geometry"] == "circle":
        electrodePositions[i] = [parameters["radius"]*np.cos(electrodes[i][0]/360*2*np.pi),parameters["radius"]*np.sin(electrodes[i][0]/360*2*np.pi)]
    


import matplotlib.pylab as plt


fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))



for i in range(len(electrodes)):
    if   i == parameters["outputElectrode"]: color= "blue"
    elif i == parameters["inputElectrode1"]: color= "red"
    elif i == parameters["inputElectrode2"]: color= "red"
    else:                                    color= "green"
    
    if parameters["geometry"] == "rect":
        ax.scatter(*electrodePositions[i],c=color,marker=".",s=400,zorder = -1)
    elif parameters["geometry"] == "circle":
        width = 8
        electrodeWidth = parameters["electrodeWidth"]/(parameters["radius"]*2*np.pi)*360 #in degrees
        ax.add_artist(Wedge((0,0),parameters["radius"]+width/2,electrodes[i][0]-electrodeWidth/2,electrodes[i][0]+electrodeWidth/2, width = width, fc=color,ec=color,zorder = -1))






# for i in range(acceptorPos.shape[0]):
#     ax.text(acceptorPos[i,0],acceptorPos[i,1], f"{i}",
#             horizontalalignment='center',
#             verticalalignment='center',
#             fontsize=10, color='red',
#             transform=ax.transData)
#     ax.add_artist(plt.Circle((acceptorPos[i,0],acceptorPos[i,1]),8,fc='none',ec="r"))

# for i in range(donorPos.shape[0]):
#     ax.text(donorPos[i,0],donorPos[i,1], f"{i}",
#             horizontalalignment='center',
#             verticalalignment='center',
#             fontsize=10, color='blue',
#             transform=ax.transData)
#     ax.add_artist(plt.Circle((donorPos[i,0],donorPos[i,1]),8,fc='none',ec="b"))

ax.scatter(acceptorPos[:,0],acceptorPos[:,1],c="k",marker=".",s=20)
ax.scatter(donorPos[:,0],donorPos[:,1],c="k",marker="x",s=20)

if parameters["geometry"] == "circle":
    ax.add_artist(plt.Circle((0,0),parameters["radius"],fc='none',ec="k",zorder = -2)) 
    ax.axis('off')
    
if parameters["geometry"] == "rect":
    ax.set_xlim(0,parameters["lenX"])
    ax.set_ylim(0,parameters["lenY"])
elif parameters["geometry"] == "circle":
    ax.set_xlim(-parameters["radius"]*1.1,parameters["radius"]*1.1)
    ax.set_ylim(-parameters["radius"]*1.1,parameters["radius"]*1.1)
    
ax.set_aspect('equal')

ax.set_xticks([], [])
ax.set_yticks([], [])

plt.savefig(join(pathToSimFolder,"plotDevice.png"),bbox_inches="tight",dpi=300)    
plt.show()
plt.close()
fig=None
# print(parameters)
