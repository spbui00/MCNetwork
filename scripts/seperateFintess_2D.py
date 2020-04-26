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

bins = 50


print("################ XOR ################")



D_11_00 = currents[:,1,1]-currents[:,0,0]
D_10_01 = currents[:,1,0]-currents[:,0,1]
D_DIFF  = (currents[:,1,1]+currents[:,0,0]-currents[:,0,1]-currents[:,1,0])/2


D0_counts, D0_bins = np.histogram(D_11_00, bins = bins, density = True)
D1_counts, D1_bins = np.histogram(D_10_01, bins = bins, density = True)
Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)



D0_D1_counts, _ , _       = np.histogram2d(D_11_00, D_10_01, bins = bins, density = True)
D0_Dd_counts, _ , _       = np.histogram2d(D_11_00, D_DIFF,  bins = bins, density = True)
D1_Dd_counts, _ , _       = np.histogram2d(D_10_01, D_DIFF,  bins = bins, density = True)

linSeperated_D0_D1_counts = np.multiply( * np.meshgrid(D0_counts, D1_counts, indexing = "ij"))
linSeperated_D0_Dd_counts = np.multiply( * np.meshgrid(D0_counts, Dd_counts, indexing = "ij"))
linSeperated_D1_Dd_counts = np.multiply( * np.meshgrid(D1_counts, Dd_counts, indexing = "ij"))



D0_D1_counts              *= (D0_bins[1]-D0_bins[0]) * (D1_bins[1]-D1_bins[0])
D0_Dd_counts              *= (D0_bins[1]-D0_bins[0]) * (Dd_bins[1]-Dd_bins[0])
D1_Dd_counts              *= (D1_bins[1]-D1_bins[0]) * (Dd_bins[1]-Dd_bins[0])
linSeperated_D0_D1_counts *= (D0_bins[1]-D0_bins[0]) * (D1_bins[1]-D1_bins[0])
linSeperated_D0_Dd_counts *= (D0_bins[1]-D0_bins[0]) * (Dd_bins[1]-Dd_bins[0])
linSeperated_D1_Dd_counts *= (D1_bins[1]-D1_bins[0]) * (Dd_bins[1]-Dd_bins[0])






fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D0_D1_counts), np.min(D0_D1_counts))
vmax = max(np.max(linSeperated_D0_D1_counts), np.max(D0_D1_counts))

im1 = ax1.imshow(D0_D1_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D0_D1_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_density2D_D0_D1_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)








fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D0_Dd_counts), np.min(D0_Dd_counts))
vmax = max(np.max(linSeperated_D0_Dd_counts), np.max(D0_Dd_counts))

im1 = ax1.imshow(D0_Dd_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D0_Dd_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_density2D_D0_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)







fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D1_Dd_counts), np.min(D1_Dd_counts))
vmax = max(np.max(linSeperated_D1_Dd_counts), np.max(D1_Dd_counts))

im1 = ax1.imshow(D1_Dd_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D1_Dd_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_density2D_D1_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)





########################### fitness plots ###########################
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


D0_D1_fitness_seperated = np.sum(probd * fitness , axis = 2)
D0_Dd_fitness_seperated = np.sum(prob1 * fitness , axis = 1)
D1_Dd_fitness_seperated = np.sum(prob0 * fitness , axis = 0)


D0_D1_fitness = np.histogram2d(D_11_00, D_10_01, bins = bins, weights = controlFitness)[0]/np.histogram2d(D_11_00, D_10_01, bins = bins)[0]
D0_Dd_fitness = np.histogram2d(D_11_00, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_11_00, D_DIFF , bins = bins)[0]
D1_Dd_fitness = np.histogram2d(D_10_01, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_10_01, D_DIFF , bins = bins)[0]


fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D0_D1_fitness_seperated), np.min(D0_D1_fitness))
vmax = max(np.max(D0_D1_fitness_seperated), np.max(D0_D1_fitness))

im1 = ax1.imshow(D0_D1_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(D0_D1_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_fitness2D_D0_D1_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)



fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D0_Dd_fitness_seperated), np.min(D0_Dd_fitness))
vmax = max(np.max(D0_Dd_fitness_seperated), np.max(D0_Dd_fitness))

im1 = ax1.imshow(D0_Dd_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(D0_Dd_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_fitness2D_D0_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)




fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D1_Dd_fitness_seperated), np.min(D1_Dd_fitness))
vmax = max(np.max(D1_Dd_fitness_seperated), np.max(D1_Dd_fitness))

im1 = ax1.imshow(D1_Dd_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")
im2 = ax2.imshow(D1_Dd_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{2}^\textrm{X}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{X}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"XOR_fitness2D_D1_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)



print("################ AND ################")



D_00_01 = currents[:,0,0]-currents[:,0,1]
D_01_10 = currents[:,0,1]-currents[:,1,0]
D_DIFF  = currents[:,1,1] -(currents[:,0,0]+currents[:,0,1]+currents[:,1,0])/3


D0_counts, D0_bins = np.histogram(D_00_01, bins = bins, density = True)
D1_counts, D1_bins = np.histogram(D_01_10, bins = bins, density = True)
Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)



D0_D1_counts, _ , _       = np.histogram2d(D_00_01, D_01_10, bins = bins, density = True)
D0_Dd_counts, _ , _       = np.histogram2d(D_00_01, D_DIFF,  bins = bins, density = True)
D1_Dd_counts, _ , _       = np.histogram2d(D_01_10, D_DIFF,  bins = bins, density = True)

linSeperated_D0_D1_counts = np.multiply( * np.meshgrid(D0_counts, D1_counts, indexing = "ij"))
linSeperated_D0_Dd_counts = np.multiply( * np.meshgrid(D0_counts, Dd_counts, indexing = "ij"))
linSeperated_D1_Dd_counts = np.multiply( * np.meshgrid(D1_counts, Dd_counts, indexing = "ij"))



D0_D1_counts              *= (D0_bins[1]-D0_bins[0]) * (D1_bins[1]-D1_bins[0])
linSeperated_D0_D1_counts *= (D0_bins[1]-D0_bins[0]) * (D1_bins[1]-D1_bins[0])
D0_Dd_counts              *= (D0_bins[1]-D0_bins[0]) * (Dd_bins[1]-Dd_bins[0])
linSeperated_D0_Dd_counts *= (D0_bins[1]-D0_bins[0]) * (Dd_bins[1]-Dd_bins[0])
D1_Dd_counts              *= (D1_bins[1]-D1_bins[0]) * (Dd_bins[1]-Dd_bins[0])
linSeperated_D1_Dd_counts *= (D1_bins[1]-D1_bins[0]) * (Dd_bins[1]-Dd_bins[0])






fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D0_D1_counts), np.min(D0_D1_counts))
vmax = max(np.max(linSeperated_D0_D1_counts), np.max(D0_D1_counts))

im1 = ax1.imshow(D0_D1_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D0_D1_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_density2D_D0_D1_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)








fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D0_Dd_counts), np.min(D0_Dd_counts))
vmax = max(np.max(linSeperated_D0_Dd_counts), np.max(D0_Dd_counts))

im1 = ax1.imshow(D0_Dd_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D0_Dd_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_density2D_D0_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)







fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(linSeperated_D1_Dd_counts), np.min(D1_Dd_counts))
vmax = max(np.max(linSeperated_D1_Dd_counts), np.max(D1_Dd_counts))

im1 = ax1.imshow(D1_Dd_counts             , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")
im2 = ax2.imshow(linSeperated_D1_Dd_counts, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")

# fig.subplots_adjust(right=0.82)
# cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
# fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_density2D_D1_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)





########################### fitness plots ###########################
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


fitness = 1-(I_normed[:,:,:,0,0] + 1-I_normed[:,:,:,0,1] + 1-I_normed[:,:,:,1,0] + I_normed[:,:,:,1,1])/4 #XOR




#calc probabilities

prob0, prob1, probd = np.meshgrid((D0_bins[1:]-D0_bins[:-1]) * D0_counts,
                                  (D1_bins[1:]-D1_bins[:-1]) * D1_counts,
                                  (Dd_bins[1:]-Dd_bins[:-1]) * Dd_counts, indexing = "ij")


D0_D1_fitness_seperated = np.sum(probd * fitness , axis = 2)
D0_Dd_fitness_seperated = np.sum(prob1 * fitness , axis = 1)
D1_Dd_fitness_seperated = np.sum(prob0 * fitness , axis = 0)


D0_D1_fitness = np.histogram2d(D_00_01, D_01_10, bins = bins, weights = controlFitness)[0]/np.histogram2d(D_00_01, D_01_10, bins = bins)[0]
D0_Dd_fitness = np.histogram2d(D_00_01, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_00_01, D_DIFF , bins = bins)[0]
D1_Dd_fitness = np.histogram2d(D_01_10, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_01_10, D_DIFF , bins = bins)[0]


fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D0_D1_fitness_seperated), np.min(D0_D1_fitness))
vmax = max(np.max(D0_D1_fitness_seperated), np.max(D0_D1_fitness))

im1 = ax1.imshow(D0_D1_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(D0_D1_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ D1_bins[0], D1_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_fitness2D_D0_D1_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)



fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D0_Dd_fitness_seperated), np.min(D0_Dd_fitness))
vmax = max(np.max(D0_Dd_fitness_seperated), np.max(D0_Dd_fitness))

im1 = ax1.imshow(D0_Dd_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")
im2 = ax2.imshow(D0_Dd_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D0_bins[-1], D0_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{1}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_fitness2D_D0_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)




fig, (ax1, ax2)=plt.subplots(1,2,figsize=(4.980614173228346,3.2), sharey=True)

vmin = min(np.min(D1_Dd_fitness_seperated), np.min(D1_Dd_fitness))
vmax = max(np.max(D1_Dd_fitness_seperated), np.max(D1_Dd_fitness))

im1 = ax1.imshow(D1_Dd_fitness          , cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")
im2 = ax2.imshow(D1_Dd_fitness_seperated, cmap = "gnuplot", vmin = vmin, vmax = vmax, extent = [ Dd_bins[0], Dd_bins[-1], D1_bins[-1], D1_bins[0]], aspect = "auto")

ax1.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")
ax1.set_ylabel(r"$\scriptsize \Delta_{2}^\textrm{A}$")
ax2.set_xlabel(r"$\scriptsize \Delta_{3}^\textrm{A}$")

fig.subplots_adjust(right=0.82)
cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
fig.colorbar(im1, cax=cbar_ax)
plt.savefig(join(pathToSimFolder,f"AND_fitness2D_D1_Dd_bins_{bins}.png"),bbox_inches="tight",dpi=300)    
# plt.show()
plt.close(fig)


