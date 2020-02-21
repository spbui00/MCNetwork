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
    optEnergy=np.array(dataFile["/optEnergy"][:])
    try:
        generations=np.array(dataFile["/generation"][:])
        mode="genetic"
    except KeyError:
        accepted=np.array(dataFile["/accepted"][:],dtype=bool)
        notAccepted=np.invert(accepted)
        mode="MC"
        

# print(accepted)



if mode=="genetic":
    genBest=np.empty(optEnergy.shape)
    for i in range(int(optEnergy.shape[0]/25)):
        genBest[i*25:(i+1)*25]=max(optEnergy[i*25:(i+1)*25])
        

fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))



if mode=="genetic":
    ax.plot(optEnergy,".",ms=1,color="darkgreen",label="all")
    ax.plot(genBest,color="darkblue",label="gen best")

if mode=="MC":
    ax.plot(np.arange(optEnergy.shape[0])[notAccepted[:,0]],optEnergy[notAccepted],".",ms=1,color="darkred",label="not accepted")
    ax.plot(np.arange(optEnergy.shape[0])[accepted[:,0]],optEnergy[accepted],".",ms=1,color="darkgreen",label="accepted")
    


ax.plot(np.maximum.accumulate(optEnergy),color="darkorange",label="best")

# ax.set_xlim(-0.15,0.65)
ax.set_ylim(0.25,1.1)

ax.set_xlabel("iteration")
ax.set_ylabel("optEnergy")


ax.legend()


plt.savefig(join(pathToSimFolder,"convergence.png"),bbox_inches="tight",dpi=300)    
plt.show()
plt.close()
fig=None


