#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np

from time import sleep

if len(argv)>1:
    pathToSimFolder=argv[1]
else:
    pathToSimFolder="../data/"

parameters,electrodes=readParameters(pathToSimFolder)

fileOpenTries = 0
while fileOpenTries < 50:
    fileOpenTries += 1
    try:
        with h5py.File(join(pathToSimFolder,"data.hdf5"),"r") as dataFile:
            optEnergy=np.array(dataFile["/optEnergy"][:])
            try:
                generations=np.array(dataFile["/generation"][:])
                mode="genetic"
            except KeyError:
                mode="MC"
                try:
                    accepted=np.array(dataFile["/accepted"][:],dtype=bool)
                    notAccepted=np.invert(accepted)
                except KeyError:
                    accepted=np.ones(optEnergy.shape,dtype=bool)
                    notAccepted=np.invert(accepted)
        break
    except OSError as e:
        if "No such file" in repr(e) :
            raise e
        else:
            print(f"could not open file. try number {fileOpenTries}")
            sleep(1)
            
        

# print(accepted)



if mode=="genetic":
    genBest=np.empty(optEnergy.shape)
    for i in range(int(optEnergy.shape[0]/25)):
        genBest[i*25:(i+1)*25]=max(optEnergy[i*25:(i+1)*25])
        

fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


ax.plot(np.maximum.accumulate(optEnergy),color="darkorange",label="best")

if mode=="genetic":
    ax.plot(optEnergy,".",ms=1,color="darkgreen",label="all")
    ax.plot(genBest,color="darkblue",label="gen best")

if mode=="MC":
    ax.plot(np.arange(optEnergy.shape[0])[notAccepted[:,0]],optEnergy[notAccepted],".",ms=1,color="darkred",label="not accepted")
    ax.plot(np.arange(optEnergy.shape[0])[accepted[:,0]],optEnergy[accepted],".",ms=1,color="darkgreen",label="accepted")
    


# ax.set_xlim(-0.15,0.65)
ax.set_ylim(0.15,1.05)

ax.set_xlabel("iteration")
ax.set_ylabel("optEnergy")


ax.legend()


plt.savefig(join(pathToSimFolder,"convergence.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close()
fig=None


