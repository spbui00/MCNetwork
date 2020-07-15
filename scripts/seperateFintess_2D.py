#!/usr/bin/python3
from tools import *
from sys import argv
from os.path import join

import h5py
import matplotlib.pylab as plt
import numpy as np


###################################################### see old version at the end (easier to understand, but less generalized and not using pca) ########################################################################

if len(argv) > 1:
    pathToSimFolder = argv[1]
else:
    pathToSimFolder = "../data/"


# controlFitness = np.load(join(pathToSimFolder,"fitness.npy"       ))
currents = np.load(join(pathToSimFolder, "currents.npy"))
currentsUncert = np.load(join(pathToSimFolder, "currentsUncert.npy"))

min_ = np.min(currents, axis=(1, 2))
max_ = np.max(currents, axis=(1, 2))
normedCurrents = ((currents.T - min_) / (max_ - min_)).T


samples = currents.shape[0]
print("samples: ", samples)

# minFitness = 0.8
bins = 200
# fitnessBins = 100


D = np.array(
    [
        currents[:, 0, 0] - currents[:, 0, 1],
        currents[:, 0, 0] - currents[:, 1, 0],
        currents[:, 0, 0] - currents[:, 1, 1],
    ]
)
N = D.shape[0]


C = np.cov(D)
_, M = np.linalg.eig(C)

# M = np.array([[1,-1,0],[0,0,-1],[0.5,0.5,-0.5]]).T #<<<<<-------- define transformation matrix here and ignore PCA
Delta = M.T @ D

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Delta_bins, Delta_counts = [], []

for i in range(N):
    Delta_counts_, Delta_bins_ = np.histogram(Delta[i], bins=bins, density=True)
    Delta_bins.append(Delta_bins_)
    Delta_counts.append(Delta_counts_)


Delta_mesh = np.array(
    np.meshgrid(
        *[(Delta_bins[i][1:] + Delta_bins[i][:-1]) / 2 for i in range(N)], indexing="ij"
    )
)
probabilities = np.array(
    np.meshgrid(
        *[(Delta_bins[i][1:] - Delta_bins[i][:-1]) * Delta_counts[i] for i in range(N)],
        indexing="ij",
    )
)
totProbab = np.prod(probabilities, axis=0)


D_mesh = np.einsum("ji...,i...->j...", np.linalg.inv(M.T), Delta_mesh)

diffs = np.zeros((bins, bins, bins, 2, 2, 2, 2))

diffs[:, :, :, 0, 0, 0, 1] = D_mesh[0]
diffs[:, :, :, 0, 0, 1, 0] = D_mesh[1]
diffs[:, :, :, 0, 0, 1, 1] = D_mesh[2]
diffs[:, :, :, 0, 1, 1, 0] = D_mesh[1] - D_mesh[0]
diffs[:, :, :, 0, 1, 1, 1] = D_mesh[2] - D_mesh[0]
diffs[:, :, :, 1, 0, 1, 1] = D_mesh[2] - D_mesh[1]

diffs[:, :, :, 0, 1, 0, 0] = -diffs[:, :, :, 0, 0, 0, 1]
diffs[:, :, :, 1, 0, 0, 0] = -diffs[:, :, :, 0, 0, 1, 0]
diffs[:, :, :, 1, 1, 0, 0] = -diffs[:, :, :, 0, 0, 1, 1]
diffs[:, :, :, 1, 0, 0, 1] = -diffs[:, :, :, 0, 1, 1, 0]
diffs[:, :, :, 1, 1, 0, 1] = -diffs[:, :, :, 0, 1, 1, 1]
diffs[:, :, :, 1, 1, 1, 0] = -diffs[:, :, :, 1, 0, 1, 1]

max_min = np.max(np.abs(diffs), axis=(3, 4, 5, 6))

I_normed = np.zeros((bins, bins, bins, 2, 2))

I_normed[:, :, :, 0, 0] = np.max(diffs[:, :, :, 0, 0, :, :], axis=(3, 4)) / max_min
I_normed[:, :, :, 0, 1] = np.max(diffs[:, :, :, 0, 1, :, :], axis=(3, 4)) / max_min
I_normed[:, :, :, 1, 0] = np.max(diffs[:, :, :, 1, 0, :, :], axis=(3, 4)) / max_min
I_normed[:, :, :, 1, 1] = np.max(diffs[:, :, :, 1, 1, :, :], axis=(3, 4)) / max_min


def getFitness(gate, normed_currents):
    def gateFunc(in1, in2):
        if gate == "AND":
            return in1 & in2
        if gate == "NAND":
            return not (in1 & in2)
        if gate == "OR":
            return in1 | in2
        if gate == "NOR":
            return not (in1 | in2)
        if gate == "XOR":
            return in1 ^ in2
        if gate == "XNOR":
            return not (in1 ^ in2)

    if len(normed_currents.shape) == 5:
        return 1 - 0.25 * np.sum(
            [
                abs(
                    gateFunc(int(i / 2), i % 2)
                    - normed_currents[:, :, :, int(i / 2), i % 2]
                )
                for i in range(4)
            ],
            axis=0,
        )
    elif len(normed_currents.shape) == 3:
        return 1 - 0.25 * np.sum(
            [
                abs(gateFunc(int(i / 2), i % 2) - normed_currents[:, int(i / 2), i % 2])
                for i in range(4)
            ],
            axis=0,
        )


################################################ 2D distr ################################################


axlims = [[-0.01, 0.07], [-0.015, 0.015], [-0.01, 0.01]]

permutations = [(0, 1), (0, 2), (1, 2)]

counts2D = []
for i in range(len(permutations)):
    counts2D.append(
        np.histogram2d(
            Delta[permutations[i][0]],
            Delta[permutations[i][1]],
            bins=bins,
            density=True,
        )[0]
    )

counts2D_linSeperated = []
for i in range(len(permutations)):
    counts2D_linSeperated.append(
        np.multiply(
            *np.meshgrid(
                Delta_counts[permutations[i][0]],
                Delta_counts[permutations[i][1]],
                indexing="ij",
            )
        )
    )

for i in range(len(permutations)):
    counts2D[i] *= (
        Delta_bins[permutations[i][0]][1] - Delta_bins[permutations[i][0]][0]
    ) * (Delta_bins[permutations[i][1]][1] - Delta_bins[permutations[i][1]][0])
    counts2D_linSeperated[i] *= (
        Delta_bins[permutations[i][0]][1] - Delta_bins[permutations[i][0]][0]
    ) * (Delta_bins[permutations[i][1]][1] - Delta_bins[permutations[i][1]][0])


for i in range(len(permutations)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(4.980614173228346, 3.2), sharey=True)

    vmin = min(np.min(counts2D_linSeperated[i]), np.min(counts2D[i]))
    vmax = max(np.max(counts2D_linSeperated[i]), np.max(counts2D[i]))

    im1 = ax1.imshow(
        counts2D[i],
        cmap="gnuplot",
        vmin=vmin,
        vmax=vmax,
        extent=[
            Delta_bins[permutations[i][1]][0],
            Delta_bins[permutations[i][1]][-1],
            Delta_bins[permutations[i][0]][-1],
            Delta_bins[permutations[i][0]][0],
        ],
        aspect="auto",
    )
    im2 = ax2.imshow(
        counts2D_linSeperated[i],
        cmap="gnuplot",
        vmin=vmin,
        vmax=vmax,
        extent=[
            Delta_bins[permutations[i][1]][0],
            Delta_bins[permutations[i][1]][-1],
            Delta_bins[permutations[i][0]][-1],
            Delta_bins[permutations[i][0]][0],
        ],
        aspect="auto",
    )

    ax1.set_xlabel(rf"$\scriptsize \Delta_{{ {permutations[i][1]+1} }}$")
    ax1.set_ylabel(rf"$\scriptsize \Delta_{{ {permutations[i][0]+1} }}$")
    ax2.set_xlabel(rf"$\scriptsize \Delta_{{ {permutations[i][1]+1} }}$")

    ax1.set_xlim(axlims[permutations[i][1]])
    ax2.set_xlim(axlims[permutations[i][1]])
    ax1.set_ylim(axlims[permutations[i][0]])

    # fig.subplots_adjust(right=0.82)
    # cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
    # fig.colorbar(im1, cax=cbar_ax)
    plt.savefig(
        join(
            pathToSimFolder,
            f"density2D_Delta{permutations[i][0]+1}_Delta{permutations[i][1]+1}_bins_{bins}.png",
        ),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)


########################### fitness plots ###########################

gates = ["AND", "NAND", "OR", "NOR", "XOR", "XNOR"]
# gates = ["XOR"]
for gate in gates:
    fitness = getFitness(gate, I_normed)
    controlFitness = getFitness(gate, normedCurrents)

    for i in range(len(permutations)):
        fig, (ax1, ax2) = plt.subplots(
            1, 2, figsize=(4.980614173228346, 3.2), sharey=True
        )

        thirdIndex = [j for j in range(3) if j not in permutations[i]][0]
        fitness_seperated = np.sum(probabilities[thirdIndex] * fitness, axis=thirdIndex)
        fitness2D = (
            np.histogram2d(
                Delta[permutations[i][0]],
                Delta[permutations[i][1]],
                bins=bins,
                weights=controlFitness,
            )[0]
            / np.histogram2d(
                Delta[permutations[i][0]], Delta[permutations[i][1]], bins=bins
            )[0]
        )

        vmin = min(np.min(fitness_seperated), np.min(fitness2D))
        vmax = max(np.max(fitness_seperated), np.max(fitness2D))

        im1 = ax1.imshow(
            fitness2D,
            cmap="gnuplot",
            vmin=vmin,
            vmax=vmax,
            extent=[
                Delta_bins[permutations[i][1]][0],
                Delta_bins[permutations[i][1]][-1],
                Delta_bins[permutations[i][0]][-1],
                Delta_bins[permutations[i][0]][0],
            ],
            aspect="auto",
        )
        im2 = ax2.imshow(
            fitness_seperated,
            cmap="gnuplot",
            vmin=vmin,
            vmax=vmax,
            extent=[
                Delta_bins[permutations[i][1]][0],
                Delta_bins[permutations[i][1]][-1],
                Delta_bins[permutations[i][0]][-1],
                Delta_bins[permutations[i][0]][0],
            ],
            aspect="auto",
        )

        ax1.set_xlabel(rf"$\scriptsize \Delta_{{ {permutations[i][1]+1} }}$")
        ax1.set_ylabel(rf"$\scriptsize \Delta_{{ {permutations[i][0]+1} }}$")
        ax2.set_xlabel(rf"$\scriptsize \Delta_{{ {permutations[i][1]+1} }}$")

        ax1.set_xlim(axlims[permutations[i][1]])
        ax2.set_xlim(axlims[permutations[i][1]])
        ax1.set_ylim(axlims[permutations[i][0]])

        fig.subplots_adjust(right=0.82)
        cbar_ax = fig.add_axes([0.86, 0.15, 0.04, 0.7])
        fig.colorbar(im1, cax=cbar_ax)
        plt.savefig(
            join(
                pathToSimFolder,
                f"{gate}_fitness2D_Delta{permutations[i][0]+1}_Delta{permutations[i][1]+1}_bins_{bins}.png",
            ),
            bbox_inches="tight",
            dpi=300,
        )
        # plt.show()
        plt.close(fig)


######################################################  old version starts here ########################################################################
"""
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



D_10_01 = currents[:,1,0]-currents[:,0,1]
D_11_00 = currents[:,1,1]-currents[:,0,0]
D_DIFF  = (currents[:,1,1]+currents[:,0,0]-currents[:,0,1]-currents[:,1,0])/2


D0_counts, D0_bins = np.histogram(D_10_01, bins = bins, density = True)
D1_counts, D1_bins = np.histogram(D_11_00, bins = bins, density = True)
Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)



D0_D1_counts, _ , _       = np.histogram2d(D_10_01, D_11_00, bins = bins, density = True)
D0_Dd_counts, _ , _       = np.histogram2d(D_10_01, D_DIFF,  bins = bins, density = True)
D1_Dd_counts, _ , _       = np.histogram2d(D_11_00, D_DIFF,  bins = bins, density = True)

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


D0_D1_fitness = np.histogram2d(D_10_01, D_11_00, bins = bins, weights = controlFitness)[0]/np.histogram2d(D_10_01, D_11_00, bins = bins)[0]
D0_Dd_fitness = np.histogram2d(D_10_01, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_10_01, D_DIFF , bins = bins)[0]
D1_Dd_fitness = np.histogram2d(D_11_00, D_DIFF , bins = bins, weights = controlFitness)[0]/np.histogram2d(D_11_00, D_DIFF , bins = bins)[0]


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


"""
