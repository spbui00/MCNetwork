#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)
electrodeNumber=len(electrodes)

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

# print(acceptorPos)
# print(donorPos)

electrodePositions=np.empty((len(electrodes),2))

for i in range(len(electrodes)):
    if electrodes[i][1]==0: electrodePositions[i] = [0                              ,electrodes[i][0]*parameters["lenY"]]
    if electrodes[i][1]==1: electrodePositions[i] = [parameters["lenX"]             ,electrodes[i][0]*parameters["lenY"]]
    if electrodes[i][1]==2: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],0                              ]
    if electrodes[i][1]==3: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],parameters["lenY"]             ]

# print(electrodePositions)


import matplotlib.pylab as plt


fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

for i in range(len(electrodes)):
    if   i == parameters["outputElectrode"]: color= "blue"
    elif i == parameters["inputElectrode1"]: color= "red"
    elif i == parameters["inputElectrode2"]: color= "red"
    else:                                    color= "green"
    
    ax.scatter(*electrodePositions[i],c=color,marker=".",s=200)



ax.scatter(acceptorPos[:,0],acceptorPos[:,1],c="k",marker=".",s=20)
ax.scatter(donorPos[:,0],donorPos[:,1],c="k",marker="x",s=20)


ax.set_xlim(0,parameters["lenX"])
ax.set_ylim(0,parameters["lenY"])

ax.set_aspect('equal')


fileNumber=4
data=np.genfromtxt(join(pathToSimFolder,f"swapTrackFile{fileNumber}.txt"),delimiter=";",dtype=int)

maxIndex=np.max(data)
added=(maxIndex+1)*data[:,0]+data[:,1]



bins=np.bincount(added)
bins.resize(maxIndex+1,maxIndex+1)


#log
# bins[np.where(bins==0)]=1
# bins=np.log(bins)

#sqrt
bins=np.sqrt(bins)

absBins=bins+bins.T
absBins=absBins/np.max(absBins)



def colorMaker(x):
    from matplotlib import colors
    from scipy.interpolate import interp1d
    
    cols=["darkred","k","darkblue"]
    
    rgbaData=np.array([colors.to_rgba(c) for c in cols])
    rInterpolater=interp1d(np.linspace(0,1,len(cols)),rgbaData[:,0])
    gInterpolater=interp1d(np.linspace(0,1,len(cols)),rgbaData[:,1])
    bInterpolater=interp1d(np.linspace(0,1,len(cols)),rgbaData[:,2])
    return np.array([rInterpolater(x),gInterpolater(x),bInterpolater(x),1])
    

# fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


# for i in range(100):
#     ax.plot([i],[1],".",color=colorMaker(i/100))


# plt.show()
# plt.close()
# fig=None







for i in range(bins.shape[0]):
    if i>=parameters["acceptorNumber"]:
        x1,y1=electrodePositions[i-int(parameters["acceptorNumber"])][0],electrodePositions[i-int(parameters["acceptorNumber"])][1]
    else:
        x1,y1=acceptorPos[i,0],acceptorPos[i,1]
        
    for j in range(i):
        if j>=parameters["acceptorNumber"]:
            x2,y2=electrodePositions[j-int(parameters["acceptorNumber"])][0],electrodePositions[j-int(parameters["acceptorNumber"])][1]
        else:
            x2,y2=acceptorPos[j,0],acceptorPos[j,1]

        #which index is next to center?
        if np.sqrt(x1**2+y1**2)>np.sqrt(x2**2+y2**2):
            index1=i
            index2=j
        else:
            index1=j
            index2=i

        # ax.plot([x1,x2],[y1,y2],"k-",alpha=bins[i,j])
        if (bins[i,j]+bins[j,i]) !=0:
            ax.plot([x1,x2],[y1,y2],"-",alpha=absBins[i,j],color=colorMaker(bins[index1,index2]/(bins[i,j]+bins[j,i])),linewidth=5)

ax.set_xlim(-0.01*parameters["lenX"],1.01*parameters["lenX"])
ax.set_ylim(-0.01*parameters["lenY"],1.01*parameters["lenY"])


plt.savefig(join(pathToSimFolder,f"swapTrack{fileNumber}.png"),bbox_inches="tight",dpi=300)    

# plt.show()
plt.close()
fig=None





# acceptorPos=np.zeros((int(parameters["acceptorNumber"]),2))
# donorPos   =np.zeros((int(parameters["acceptorNumber"]*parameters["compensationFactor"]),2))
# with open(join(pathToSimFolder,"device.txt")) as deviceFile:
#     line=next(deviceFile)
#     for i in range(acceptorPos.shape[0]):
#         acceptorPos[i]=next(deviceFile).split(" ")
#     line=next(deviceFile)
#     line=next(deviceFile)
#     for i in range(donorPos.shape[0]):
#         donorPos[i]=next(deviceFile).split(" ")

# # print(acceptorPos)
# # print(donorPos)

# electrodePositions=np.empty((len(electrodes),2))

# for i in range(len(electrodes)):
#     if electrodes[i][1]==0: electrodePositions[i] = [0                              ,electrodes[i][0]*parameters["lenY"]]
#     if electrodes[i][1]==1: electrodePositions[i] = [parameters["lenX"]             ,electrodes[i][0]*parameters["lenY"]]
#     if electrodes[i][1]==2: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],0                              ]
#     if electrodes[i][1]==3: electrodePositions[i] = [electrodes[i][0]*parameters["lenX"],parameters["lenY"]             ]

# # print(electrodePositions)


# import matplotlib.pylab as plt


# fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

# for i in range(len(electrodes)):
#     if   i == parameters["outputElectrode"]: color= "blue"
#     elif i == parameters["inputElectrode1"]: color= "red"
#     elif i == parameters["inputElectrode2"]: color= "red"
#     else:                                    color= "green"
    
#     ax.scatter(*electrodePositions[i],c=color,marker=".",s=200)



# ax.scatter(acceptorPos[:,0],acceptorPos[:,1],c="k",marker=".",s=20)
# ax.scatter(donorPos[:,0],donorPos[:,1],c="k",marker="x",s=20)


# ax.set_xlim(0,parameters["lenX"])
# ax.set_ylim(0,parameters["lenY"])

# ax.set_aspect('equal')


# fileNumber=4
# data=np.genfromtxt(join(pathToSimFolder,f"swapTrackFile{fileNumber}.txt"),delimiter=";",dtype=int)

# maxIndex=np.max(data)
# added=(maxIndex+1)*data[:,0]+data[:,1]



# bins=np.bincount(added)
# bins.resize(maxIndex+1,maxIndex+1)


# #log
# # bins[np.where(bins==0)]=1
# # bins=np.log(bins)

# #sqrt
# bins=np.sqrt(bins)


# bins=bins/np.max(bins)


# for i in range(bins.shape[0]):
#     if i>=parameters["acceptorNumber"]:
#         x1,y1=electrodePositions[i-int(parameters["acceptorNumber"])][0],electrodePositions[i-int(parameters["acceptorNumber"])][1]
#     else:
#         x1,y1=acceptorPos[i,0],acceptorPos[i,1]
        
#     for j in range(bins.shape[0]):
#         if j>=parameters["acceptorNumber"]:
#             x2,y2=electrodePositions[j-int(parameters["acceptorNumber"])][0],electrodePositions[j-int(parameters["acceptorNumber"])][1]
#         else:
#             x2,y2=acceptorPos[j,0],acceptorPos[j,1]

#         ax.plot([x1,x2],[y1,y2],"k-",alpha=bins[i,j])

# ax.set_xlim(-0.01*parameters["lenX"],1.01*parameters["lenX"])
# ax.set_ylim(-0.01*parameters["lenY"],1.01*parameters["lenY"])


# plt.savefig(join(pathToSimFolder,f"swapTrack{fileNumber}.png"),bbox_inches="tight",dpi=300)    

# # plt.show()
# plt.close()
# fig=None













# data=np.genfromtxt(join(pathToSimFolder,"swapTrackFile1.txt"),delimiter=";",dtype=int)

# maxIndex=np.max(data)
# added=(maxIndex+1)*data[:,0]+data[:,1]


# fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

# bins=np.bincount(added)

# # for i in range(len(bins)):
#     # print(i%,bins[i])

# bins.resize(maxIndex+1,maxIndex+1)
# im=ax.imshow(bins,cmap="Greys")

# # ax.set_xticks(np.arange(0,parameters["acceptorNumber"]+electrodeNumber,2))
# # ax.set_xticklabels([f"A{i:.0f}" for i in np.arange(0,parameters["acceptorNumber"],2)]+[f"E{i}" for i in np.arange(0,electrodeNumber,2)])
# # ax.set_yticks(np.arange(0,parameters["acceptorNumber"]+electrodeNumber,2))
# # ax.set_yticklabels([f"A{i:.0f}" for i in np.arange(0,parameters["acceptorNumber"],2)]+[f"E{i}" for i in np.arange(0,electrodeNumber,2)])

# # plt.savefig(join(pathToSimFolder,"fitness_1D.png"),bbox_inches="tight",dpi=300)    

# plt.colorbar(im)

# plt.show()
# plt.close()
# fig=None

