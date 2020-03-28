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
            voltages=np.array(dataFile["/voltages"][:])
            optEnergy=np.array(dataFile["/optEnergy"][:])
            while True:
                try:
                    generations=np.array(dataFile["/generation"][:])
                    mode="genetic"
                    break
                except KeyError:
                    pass
                try:
                    basinAccepted=np.array(dataFile["/basinAccepted"][:],dtype=int)
                    accepted=basinAccepted.astype(bool)
                    notAccepted=np.invert(accepted)
                    mode="basinHop"
                    break
                except KeyError:
                    pass
                mode="MC"
                try:
                    accepted=np.array(dataFile["/accepted"][:],dtype=bool)
                    notAccepted=np.invert(accepted)
                except KeyError: 
                    accepted=np.ones(optEnergy.shape,dtype=bool) #support for deprecated version
                    notAccepted=np.invert(accepted)
                break
        break
    except OSError as e:
        if "No such file" in repr(e) :
            raise e
        else:
            print(f"could not open file. try number {fileOpenTries}")
            sleep(1)
            



cotrolElectrodeIndices = list(range(0,len(electrodes)))
cotrolElectrodeIndices.remove(parameters["outputElectrode"])
cotrolElectrodeIndices.remove(parameters["inputElectrode1"])
cotrolElectrodeIndices.remove(parameters["inputElectrode2"])



cotrolVoltages = voltages[:,cotrolElectrodeIndices]


distance = 0
meanRange = 1000

displace = []

for i in range(int(distance+meanRange/2),cotrolVoltages.shape[0]):
    mean = np.mean(cotrolVoltages[int(i-distance-meanRange/2):int(i-distance+meanRange/2), :],axis = 0)
    # displace.append(np.sqrt(np.sum((cotrolVoltages[i])**2)))
    displace.append(np.sqrt(np.sum((mean-cotrolVoltages[i])**2)))



MSD = np.sum((cotrolVoltages[0]-cotrolVoltages[:])**2,axis=1)

fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


ax.plot(range(len(MSD)), MSD, "r-",label = "MSD")
ax2=ax.twinx()
ax2.plot(range(int(distance+meanRange/2),cotrolVoltages.shape[0]), displace, "k-",label="displacement")

ax.legend()
ax2.legend()

ax.set_xlabel("step")
ax.set_ylabel("displacement")

plt.savefig(join(pathToSimFolder,"displacement.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close()
fig=None





# print(mode)



if mode=="genetic":
    genBest=np.empty(optEnergy.shape)
    for i in range(int(optEnergy.shape[0]/25)):
        genBest[i*25:(i+1)*25]=max(optEnergy[i*25:(i+1)*25])
        

fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


ax.plot(np.maximum.accumulate(optEnergy),color="darkorange",label="best")

if mode=="genetic":
    ax.plot(optEnergy,".",ms=1,color="darkgreen",label="all")
    ax.plot(genBest,color="darkblue",label="gen best")

if mode=="MC" or mode =="basinHop":
    ax.plot(np.arange(optEnergy.shape[0])[notAccepted[:,0]],optEnergy[notAccepted],".",ms=1,color="darkred",label="not accepted")
    ax.plot(np.arange(optEnergy.shape[0])[accepted[:,0]],optEnergy[accepted],".",ms=1,color="darkgreen",label="accepted")
    
if mode == "basinHop":
    for i in np.where(basinAccepted[:,0] == 2)[0]:
        ax.axvline(i, color = "darkred", zorder= -1 )
    for i in np.where(basinAccepted[:,0] == 3)[0]:
        ax.axvline(i, color = "darkgreen", zorder= -1 )
        
    basinChanges = np.where(np.logical_or(basinAccepted[:,0] == 2,  basinAccepted[:,0] == 3))[0] 
    ax.plot(np.arange(0,basinChanges[0]),np.maximum.accumulate(optEnergy[:basinChanges[0]]),color="darkblue",label="basin best")
    for i in range(1,len(basinChanges)):
        ax.plot(np.arange(basinChanges[i-1],basinChanges[i]),np.maximum.accumulate(optEnergy[basinChanges[i-1]:basinChanges[i]]),color="darkblue")
    ax.plot(np.arange(basinChanges[-1],len(optEnergy)),np.maximum.accumulate(optEnergy[basinChanges[-1]:]),color="darkblue")
        
        
        

# ax.set_xlim(-0.15,0.65)
ax.set_ylim(0.15,1.05)

ax.set_xlabel("iteration")
ax.set_ylabel("optEnergy")


ax.legend()


plt.savefig(join(pathToSimFolder,"convergence.png"),bbox_inches="tight",dpi=300)    
plt.show()
plt.close()
fig=None




fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


ax.plot(np.maximum.accumulate(optEnergy),color="darkorange",label="best")

if mode=="genetic":
    ax.plot(optEnergy,".",ms=1,color="darkgreen",label="all")
    ax.plot(genBest,color="darkblue",label="gen best")


if mode=="MC" or mode =="basinHop":
    ax.plot(np.arange(optEnergy.shape[0])[notAccepted[:,0]],optEnergy[notAccepted],".",ms=1,color="darkred",label="not accepted")
    ax.plot(np.arange(optEnergy.shape[0])[accepted[:,0]],optEnergy[accepted],".",ms=1,color="darkgreen",label="accepted")
    
if mode == "basinHop":
    for i in np.where(basinAccepted[:,0] == 2)[0]:
        ax.axvline(i, color = "darkred", zorder= -1 )
    for i in np.where(basinAccepted[:,0] == 3)[0]:
        ax.axvline(i, color = "darkgreen", zorder= -1 )
        
    basinChanges = np.where(np.logical_or(basinAccepted[:,0] == 2,  basinAccepted[:,0] == 3))[0] 
    ax.plot(np.arange(0,basinChanges[0]),np.maximum.accumulate(optEnergy[:basinChanges[0]]),color="darkblue",label="basin best")
    for i in range(1,len(basinChanges)):
        ax.plot(np.arange(basinChanges[i-1],basinChanges[i]),np.maximum.accumulate(optEnergy[basinChanges[i-1]:basinChanges[i]]),color="darkblue")
    ax.plot(np.arange(basinChanges[-1],len(optEnergy)),np.maximum.accumulate(optEnergy[basinChanges[-1]:]),color="darkblue")
        
        

ax2=ax.twinx()
ax.set_zorder(ax2.get_zorder()+1)
ax.patch.set_visible(False)

ax2.plot(range(int(distance+meanRange/2),cotrolVoltages.shape[0]), displace, "k-",label="displacement")

ax.set_ylim(0.15,1.05)

ax.set_xlabel("iteration")
ax.set_ylabel("optEnergy")
ax2.set_ylabel("displacement")


# ax.legend([line],[line.get_label()])
# ax2.legend()

plt.savefig(join(pathToSimFolder,"convergence_displacement.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close()
fig=None




