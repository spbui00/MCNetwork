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

with h5py.File(join(pathToSimFolder,"data.hdf5"),"r") as dataFile:
    currents=np.array(dataFile["/outputCurrent"][:])
    sigma=np.array(dataFile["/outputCurrentUncert"][:])
    fitness=np.array(dataFile["/fitness"][:])
    sigmaFitness=np.array(dataFile["/fitnessUncert"][:])
    voltages=np.array(dataFile["/voltages"][:])
    optEnergy=np.array(dataFile["/optEnergy"][:])

voltageScanPointsNumber = int((parameters["voltageScanMax"]-parameters["voltageScanMin"])/parameters["voltageScanResolution"]+1)
parameters["voltageScanMax"]=parameters["voltageScanMin"]+(voltageScanPointsNumber-1)*parameters["voltageScanResolution"]

currents.resize((currents.shape[0],voltageScanPointsNumber,voltageScanPointsNumber))
sigma   .resize((sigma   .shape[0],voltageScanPointsNumber,voltageScanPointsNumber))


best=np.argmax(optEnergy,axis=0)
# best=[0]
print("best at",best,"optEnergy: ",optEnergy[best],"fitness: ",fitness[best]," +-", sigmaFitness[best])
print(voltages[best])

# print(currents[best])




if voltageScanPointsNumber==2:
    fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


    ax.errorbar(np.arange(4),currents[best].flatten(),yerr=sigma[best].flatten(),marker="x",color="k",linestyle="")

    ax.set_xticks([0,1,2,3]) 
    ax.set_xticklabels([r"[0,0]",r"[1,0]",r"[0,1]",r"[1,1]"])

    plt.savefig(join(pathToSimFolder,"fitness_1D.png"),bbox_inches="tight",dpi=300)    

    plt.show()
    plt.close()
    fig=None




falseIndex=int((0  -parameters["voltageScanMin"])/parameters["voltageScanResolution"])
trueIndex =int((0.5-parameters["voltageScanMin"])/parameters["voltageScanResolution"])


maxTrue=-np.inf
minTrue=np.inf

currents2x2=np.zeros((2,2))
currentUncert2x2=np.zeros((2,2))

for b1 in [0,1]:
    for b2 in [0,1]:
        if b1: i1= trueIndex
        else : i1 =falseIndex
        if b2: i2= trueIndex
        else : i2 =falseIndex
        
        currents2x2[b1,b2]=currents[best[0],i1,i2]
        currentUncert2x2[b1,b2]=sigma[best[0],i1,i2]
        
        # print(b1,b2,logical(b1,b2,parameters["gate"]))
        if logical(b1,b2,parameters["gate"]):
            if maxTrue < currents[best[0],i1,i2]:
                maxTrue = currents[best[0],i1,i2]
        else:
            if minTrue > currents[best[0],i1,i2]:
                minTrue = currents[best[0],i1,i2]
                


fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

im=ax.imshow((currents[best[0],::-1,:]-minTrue)/(maxTrue-minTrue),cmap="Greys",extent=[parameters["voltageScanMin"]-parameters["voltageScanResolution"]/2,parameters["voltageScanMax"]+parameters["voltageScanResolution"]/2,
                                                                                     parameters["voltageScanMin"]-parameters["voltageScanResolution"]/2,parameters["voltageScanMax"]+parameters["voltageScanResolution"]/2])

ax.scatter([0,0,0.5,0.5],[0,0.5,0,0.5],c="r",marker="x")

# ax.set_xlim(-0.15,0.65)
# ax.set_ylim(-0.15,0.65)

plt.colorbar(im)

plt.savefig(join(pathToSimFolder,"fitness_2D.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close()
fig=None


if voltageScanPointsNumber!=2:
    fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


    ax.errorbar(np.arange(4),currents2x2.flatten(),yerr=currentUncert2x2.flatten(),marker="x",color="k",linestyle="")

    ax.set_xticks([0,1,2,3]) 
    ax.set_xticklabels([r"[0,0]",r"[1,0]",r"[0,1]",r"[1,1]"])

    plt.savefig(join(pathToSimFolder,"fitness_1D.png"),bbox_inches="tight",dpi=300)    

    plt.show()
    plt.close()
    fig=None

