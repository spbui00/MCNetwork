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
    

controlFitness = np.load(join(pathToSimFolder,"fitness.npy"       ))
currents       = np.load(join(pathToSimFolder,"currents.npy"      ))
currentsUncert = np.load(join(pathToSimFolder,"currentsUncert.npy"))
N = currents.shape[0]

minFitness = 0.8
relUncertThres = 10000
bins = 50


print("################ XOR ################")



D_10_01 = currents[:,1,0]-currents[:,0,1]
D_11_00 = currents[:,1,1]-currents[:,0,0]
D_DIFF  = (currents[:,1,1]+currents[:,0,0]-currents[:,0,1]-currents[:,1,0])/2

#relative uncerttainty

print("meanUncert D_10_01: ",np.mean(np.sqrt(currentsUncert[:,1,0]**2 + currentsUncert[:,0,1]**2)))
print("meanUncert D_11_00: ",np.mean(np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2)))
print("meanUncert D_DIFF: ",np.mean(0.5 * np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2 + currentsUncert[:,0,1]**2 + currentsUncert[:,1,0]**2)))

u_D_10_01 = np.sqrt(currentsUncert[:,1,0]**2 + currentsUncert[:,0,1]**2)                                                               / D_10_01
u_D_11_00 = np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2)                                                               / D_11_00
u_D_DIFF  = 0.5 * np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2 + currentsUncert[:,0,1]**2 + currentsUncert[:,1,0]**2)   / D_DIFF


valid = np.where(np.logical_and(np.logical_and(u_D_11_00 < relUncertThres, u_D_10_01 < relUncertThres), u_D_DIFF < relUncertThres) )

D_10_01 = D_10_01[valid]
D_11_00 = D_11_00[valid]
D_DIFF  = D_DIFF [valid]

n=D_11_00.shape[0]
print(f"{n} == {n/N:.2%} valid")


r = np.corrcoef(D_11_00,D_10_01)[0,1]
print(f"D_11_00 vs D_10_01: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_11_00,D_DIFF)[0,1]
print(f"D_11_00 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_10_01,D_DIFF)[0,1]
print(f"D_10_01 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
# r = np.corrcoef(currents[:,1,1],currents[:,0,0])[0,1]
# print(f"currents[:,1,1] vs currents[:,0,0]: r = {r}, T(r) = {r * np.sqrt(N-2)/np.sqrt(1-r**2)}")



D0_counts, D0_bins = np.histogram(D_10_01, bins = bins, density = True)
D1_counts, D1_bins = np.histogram(D_11_00, bins = bins, density = True)
Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)


D0, D1, Dd = np.meshgrid((D0_bins[1:]+D0_bins[:-1])/2,
                         (D1_bins[1:]+D1_bins[:-1])/2,
                         (Dd_bins[1:]+Dd_bins[:-1])/2, indexing = "ij")


diffs = np.zeros((bins,bins,bins,2,2,2,2))



diffs[:,:,:,0,0,0,1] = Dd+(+D0-D1)/2
diffs[:,:,:,0,0,1,0] = Dd+(-D0-D1)/2
diffs[:,:,:,1,1,0,1] = Dd+(+D0+D1)/2
diffs[:,:,:,1,1,1,0] = Dd+(-D0+D1)/2
diffs[:,:,:,1,0,0,1] = D0
diffs[:,:,:,1,1,0,0] = D1
diffs[:,:,:,0,1,0,0] = -diffs[:,:,:,0,0,0,1]
diffs[:,:,:,1,0,0,0] = -diffs[:,:,:,0,0,1,0]
diffs[:,:,:,0,1,1,1] = -diffs[:,:,:,1,1,0,1]
diffs[:,:,:,1,0,1,1] = -diffs[:,:,:,1,1,1,0]
diffs[:,:,:,0,1,1,0] = -diffs[:,:,:,1,0,0,1]
diffs[:,:,:,0,0,1,1] = -diffs[:,:,:,1,1,0,0]


Delta = np.max(np.abs(diffs),axis=(3,4,5,6))

I_normed = np.zeros((bins,bins,bins,2,2))

I_normed[:,:,:,0,0] = np.max(diffs[:,:,:,0,0,:,:],axis=(3,4))/Delta
I_normed[:,:,:,0,1] = np.max(diffs[:,:,:,0,1,:,:],axis=(3,4))/Delta
I_normed[:,:,:,1,0] = np.max(diffs[:,:,:,1,0,:,:],axis=(3,4))/Delta
I_normed[:,:,:,1,1] = np.max(diffs[:,:,:,1,1,:,:],axis=(3,4))/Delta


fitness = 1-(I_normed[:,:,:,0,0] + 1-I_normed[:,:,:,0,1] + 1-I_normed[:,:,:,1,0] + I_normed[:,:,:,1,1])/4 #XOR




#calc probabilities

prob0, prob1, probd = np.meshgrid((D0_bins[1:]-D0_bins[:-1]) * D0_counts,
                                  (D1_bins[1:]-D1_bins[:-1]) * D1_counts,
                                  (Dd_bins[1:]-Dd_bins[:-1]) * Dd_counts, indexing = "ij")


totProb = prob0 * prob1 * probd




indices = np.where(fitness>minFitness)
#for i in range(indices[0].shape[0]):
#    print(D0[indices[0][i],indices[1][i],indices[2][i]],D1[indices[0][i],indices[1][i],indices[2][i]],Dd[indices[0][i],indices[1][i],indices[2][i]])




print("prob:",np.sum(totProb[indices]))
print("newMean:",np.sum(totProb[indices]*fitness[indices])/np.sum(totProb[indices]))



fitnessBins = 50
fitness_counts    , fitness_bins     = np.histogram(fitness.flatten(), weights = totProb.flatten(), bins = fitnessBins, density = True)
realFitness_counts, realFitness_bins = np.histogram(controlFitness, bins = fitnessBins, density = True)


fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

ax.hist(fitness_bins[:-1]    , fitness_bins    , weights=fitness_counts    , color = color(0,2), histtype = "step", label = r"estimated")
ax.hist(realFitness_bins[:-1], realFitness_bins, weights=realFitness_counts, color = color(1,2), histtype = "step", label = r"real")

# ax.set_xlim(0.4,1)
# ax.set_ylim(0.4,1)

ax.set_xlabel(r"$\mathcal{F}$")
ax.set_ylabel(r"$P(\mathcal{F})$")

ax.legend()
plt.savefig(join(pathToSimFolder,f"fitnessDistr_XOR.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)




fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


# im = ax.imshow(np.mean(fitness, axis = 2), cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])
im = ax.imshow(fitness[:,:,2]     , cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])
# im = ax.imshow(D1[:,:,5]     , cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])

# ax.set_xlim(0.4,1)
# ax.set_ylim(0.4,1)

ax.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{X}$")

plt.colorbar(im)
plt.savefig(join(pathToSimFolder,f"mapXOR.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)





fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

ax.hist(D0_bins[:-1], D0_bins, weights=D0_counts, color = color(0,3), histtype = "step", label = r"$\scriptsize \Delta_{1}^\textrm{X}$")
ax.hist(D1_bins[:-1], D1_bins, weights=D1_counts, color = color(1,3), histtype = "step", label = r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax.hist(Dd_bins[:-1], Dd_bins, weights=Dd_counts, color = color(2,3), histtype = "step", label = r"$\scriptsize \Delta_{3}^\textrm{X}$")


np.save(join(pathToSimFolder,"XOR_D0_bins.npy"  ),D0_bins)
np.save(join(pathToSimFolder,"XOR_D0_counts.npy"),D0_counts)
np.save(join(pathToSimFolder,"XOR_D1_bins.npy"  ),D1_bins)
np.save(join(pathToSimFolder,"XOR_D1_counts.npy"),D1_counts)
np.save(join(pathToSimFolder,"XOR_Dd_bins.npy"  ),Dd_bins)
np.save(join(pathToSimFolder,"XOR_Dd_counts.npy"),Dd_counts)

ax.set_xlabel(r"$\Delta_{i}$")
ax.set_ylabel(r"$P(\Delta_{i})$")

ax.legend()

plt.savefig(join(pathToSimFolder,f"distXOR.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)

print("################ AND ################")





D_00_01 = currents[:,0,0]-currents[:,0,1]
D_01_10 = currents[:,0,1]-currents[:,1,0]
D_DIFF  = currents[:,1,1] -(currents[:,0,0]+currents[:,0,1]+currents[:,1,0])/3

#relative uncerttainty
u_D_00_01 = np.sqrt(currentsUncert[:,0,0]**2 + currentsUncert[:,0,1]**2)                                                               / D_00_01
u_D_01_10 = np.sqrt(currentsUncert[:,0,1]**2 + currentsUncert[:,1,0]**2)                                                               / D_01_10
u_D_DIFF  = np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2/9 + currentsUncert[:,0,1]**2/9 + currentsUncert[:,1,0]**2/9)   / D_DIFF


print("meanUncert D_00_01: ",np.mean(np.sqrt(currentsUncert[:,0,0]**2 + currentsUncert[:,0,1]**2)))
print("meanUncert D_01_10: ",np.mean(np.sqrt(currentsUncert[:,0,1]**2 + currentsUncert[:,1,0]**2)))
print("meanUncert D_DIFF: " ,np.mean(np.sqrt(currentsUncert[:,1,1]**2 + currentsUncert[:,0,0]**2/9 + currentsUncert[:,0,1]**2/9 + currentsUncert[:,1,0]**2/9)))


valid = np.where(np.logical_and(np.logical_and(u_D_11_00 < relUncertThres, u_D_10_01 < relUncertThres), u_D_DIFF < relUncertThres) )

D_00_01 = D_00_01[valid]
D_01_10 = D_01_10[valid]
D_DIFF  = D_DIFF [valid]


n=D_00_01.shape[0]
print(f"{n} == {n/N:.2%} valid")


r = np.corrcoef(D_00_01,D_01_10)[0,1]
print(f"D_00_01 vs D_01_10: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_00_01,D_DIFF)[0,1]
print(f"D_00_01 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_01_10,D_DIFF)[0,1]
print(f"D_01_10 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")



D0_counts, D0_bins = np.histogram(D_00_01, bins = bins, density = True)
D1_counts, D1_bins = np.histogram(D_01_10, bins = bins, density = True)
Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)


D0, D1, Dd = np.meshgrid((D0_bins[1:]+D0_bins[:-1])/2,
                         (D1_bins[1:]+D1_bins[:-1])/2,
                         (Dd_bins[1:]+Dd_bins[:-1])/2, indexing = "ij")


diffs = np.zeros((bins,bins,bins,2,2,2,2))



diffs[:,:,:,0,0,0,1] = D0
diffs[:,:,:,0,0,1,0] = D0 + D1
diffs[:,:,:,0,0,1,1] = 1/3*D1 + 2/3*D0 - Dd
diffs[:,:,:,0,1,1,0] = D1
diffs[:,:,:,0,1,1,1] = 1/3*D1 - 1/3*D0 - Dd
diffs[:,:,:,1,0,1,1] =-2/3*D1 - 1/3*D0 - Dd
diffs[:,:,:,0,1,0,0] = -diffs[:,:,:,0,0,0,1]
diffs[:,:,:,1,0,0,0] = -diffs[:,:,:,0,0,1,0]
diffs[:,:,:,1,1,0,0] = -diffs[:,:,:,0,0,1,1]
diffs[:,:,:,1,0,0,1] = -diffs[:,:,:,0,1,1,0]
diffs[:,:,:,1,1,0,1] = -diffs[:,:,:,0,1,1,1]
diffs[:,:,:,1,1,1,0] = -diffs[:,:,:,1,0,1,1]


Delta = np.max(np.abs(diffs),axis=(3,4,5,6))

I_normed = np.zeros((bins,bins,bins,2,2))

I_normed[:,:,:,0,0] = np.max(diffs[:,:,:,0,0,:,:],axis=(3,4))/Delta
I_normed[:,:,:,0,1] = np.max(diffs[:,:,:,0,1,:,:],axis=(3,4))/Delta
I_normed[:,:,:,1,0] = np.max(diffs[:,:,:,1,0,:,:],axis=(3,4))/Delta
I_normed[:,:,:,1,1] = np.max(diffs[:,:,:,1,1,:,:],axis=(3,4))/Delta


fitness = 1-(I_normed[:,:,:,0,0] + I_normed[:,:,:,0,1] + I_normed[:,:,:,1,0] + 1-I_normed[:,:,:,1,1])/4 #XOR




#calc probabilities

prob0, prob1, probd = np.meshgrid((D0_bins[1:]-D0_bins[:-1]) * D0_counts,
                                  (D1_bins[1:]-D1_bins[:-1]) * D1_counts,
                                  (Dd_bins[1:]-Dd_bins[:-1]) * Dd_counts, indexing = "ij")


totProb = prob0 * prob1 * probd




indices = np.where(fitness>minFitness)
#for i in range(indices[0].shape[0]):
#    print(D0[indices[0][i],indices[1][i],indices[2][i]],D1[indices[0][i],indices[1][i],indices[2][i]],Dd[indices[0][i],indices[1][i],indices[2][i]])




print("prob:",np.sum(totProb[indices]))
print("newMean:",np.sum(totProb[indices]*fitness[indices])/np.sum(totProb[indices]))



min_ = np.min(currents, axis = (1,2))
max_ = np.max(currents, axis = (1,2))
controlFitnessAND = 1 - 0.25 * ((currents[:,0,0]-min_)/(max_-min_) + (currents[:,1,0]-min_)/(max_-min_) + (currents[:,0,1]-min_)/(max_-min_) + 1 - (currents[:,1,1]-min_)/(max_-min_)   )



fitnessBins = 50
fitness_counts    , fitness_bins     = np.histogram(fitness.flatten(), weights = totProb.flatten(), bins = fitnessBins, density = True)
realFitness_counts, realFitness_bins = np.histogram(controlFitnessAND, bins = fitnessBins, density = True)


fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

ax.hist(fitness_bins[:-1]    , fitness_bins    , weights=fitness_counts    , color = color(0,2), histtype = "step", label = r"estimated")
ax.hist(realFitness_bins[:-1], realFitness_bins, weights=realFitness_counts, color = color(1,2), histtype = "step", label = r"real")

# ax.set_xlim(0.4,1)
# ax.set_ylim(0.4,1)

ax.set_xlabel(r"$\mathcal{F}$")
ax.set_ylabel(r"$P(\mathcal{F})$")

ax.legend()

plt.savefig(join(pathToSimFolder,f"fitnessDistr_AND.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)



fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))


# im = ax.imshow(np.mean(fitness, axis = 2), cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])
im = ax.imshow(fitness[:,:,49]     , cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])
# im = ax.imshow(D1[:,:,5]     , cmap = "gnuplot", extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]])

# ax.set_xlim(0.4,1)
# ax.set_ylim(0.4,1)

ax.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{A}$")

plt.colorbar(im)
plt.savefig(join(pathToSimFolder,f"mapAND.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)





fig, ax=plt.subplots(1,1,figsize=(4.980614173228346,3.2))

ax.hist(D0_bins[:-1], D0_bins, weights=D0_counts, color = color(0,3), histtype = "step", label = r"$\scriptsize \Delta_{1}^\textrm{A}$")
ax.hist(D1_bins[:-1], D1_bins, weights=D1_counts, color = color(1,3), histtype = "step", label = r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax.hist(Dd_bins[:-1], Dd_bins, weights=Dd_counts, color = color(2,3), histtype = "step", label = r"$\scriptsize \Delta_{3}^\textrm{A}$")

np.save(join(pathToSimFolder,"AND_D0_bins.npy"  ),D0_bins)
np.save(join(pathToSimFolder,"AND_D0_counts.npy"),D0_counts)
np.save(join(pathToSimFolder,"AND_D1_bins.npy"  ),D1_bins)
np.save(join(pathToSimFolder,"AND_D1_counts.npy"),D1_counts)
np.save(join(pathToSimFolder,"AND_Dd_bins.npy"  ),Dd_bins)
np.save(join(pathToSimFolder,"AND_Dd_counts.npy"),Dd_counts)


ax.set_xlabel(r"$\Delta_{i}$")
ax.set_ylabel(r"$P(\Delta_{i})$")


ax.legend()

plt.savefig(join(pathToSimFolder,f"distAND.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)
