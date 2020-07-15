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


currents = np.load(join(pathToSimFolder, "currents.npy"))
currentsUncert = np.load(join(pathToSimFolder, "currentsUncert.npy"))

min_ = np.min(currents, axis=(1, 2))
max_ = np.max(currents, axis=(1, 2))
normedCurrents = ((currents.T - min_) / (max_ - min_)).T


samples = currents.shape[0]
print("samples: ", samples)

minFitness = 0.8
bins = 200
fitnessBins = 100


D = np.array(
    [
        currents[:, 0, 0] - currents[:, 0, 1],
        currents[:, 0, 0] - currents[:, 1, 0],
        currents[:, 0, 0] - currents[:, 1, 1],
    ]
)
N = D.shape[0]


C = np.cov(D)
M = np.linalg.eig(C)[1].T

# M = np.array([[1,-1,0],[0,0,-1],[0.5,0.5,-0.5]]) #<<<<<-------- define transformation matrix here and ignore PCA

###print new coordinates
for i in range(3):
    print(
        f"D{i+1} = {M[i,0]+M[i,1]+M[i,2]:5.2f} I_00 +  {-M[i,0]:5.2f} I_01 +  {-M[i,1]:5.2f} I_10 +  {-M[i,2]:5.2f} I_11"
    )


print("transformation matrix:\n", M)
Delta = M @ D
print("corrcoef matrix:\n", np.corrcoef(Delta))
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Delta_bins, Delta_counts = [], []

for i in range(N):
    Delta_counts_, Delta_bins_ = np.histogram(Delta[i], bins=bins, density=True)
    Delta_bins.append(Delta_bins_)
    Delta_counts.append(Delta_counts_)

    np.save(join(pathToSimFolder, f"Delta{i}_bins.npy"), Delta_bins_)
    np.save(join(pathToSimFolder, f"Delta{i}_counts.npy"), Delta_counts_)

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


D_mesh = np.einsum("ji...,i...->j...", np.linalg.inv(M), Delta_mesh)

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


############################### delta distr ###############################

fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))


for i in range(N):
    ax.hist(
        Delta_bins[i][:-1],
        Delta_bins[i],
        weights=Delta_counts[i],
        color=color(i, N),
        histtype="step",
        label=rf"$\scriptsize \Delta_{{ {i+1} }}$",
    )

ax.set_xlabel(r"$\Delta_{i}$")
ax.set_ylabel(r"$P(\Delta_{i})$")

ax.legend()
ax.set_xlim(-0.02, 0.065)

plt.savefig(join(pathToSimFolder, f"deltaDistr.png"), bbox_inches="tight", dpi=300)
# plt.show()
plt.close(fig)

############################### fitness distr ###############################
gates = ["AND", "NAND", "OR", "NOR", "XOR", "XNOR"]
# gates = ["XOR"]
for gate in gates:
    estimatedFitness = getFitness(gate, I_normed)
    realFitness = getFitness(gate, normedCurrents)
    hitIndices = np.where(estimatedFitness > minFitness)
    hitIndices = np.where(estimatedFitness > minFitness)
    print(
        f"{gate}: {np.sum(totProbab[hitIndices]):%}    {np.sum(realFitness>minFitness)/samples:%}"
    )

    fitness_counts, fitness_bins = np.histogram(
        estimatedFitness.flatten(),
        weights=totProbab.flatten(),
        bins=fitnessBins,
        density=True,
    )
    realFitness_counts, realFitness_bins = np.histogram(
        realFitness, bins=fitnessBins, density=True
    )

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        fitness_bins[:-1],
        fitness_bins,
        weights=fitness_counts,
        color=color(0, 2),
        histtype="step",
        label=r"estimated",
    )
    ax.hist(
        realFitness_bins[:-1],
        realFitness_bins,
        weights=realFitness_counts,
        color=color(1, 2),
        histtype="step",
        label=r"real",
    )

    # ax.set_xlim(0.4,1)
    # ax.set_ylim(0,1)

    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$P(f)$")

    ax.legend(loc="best")

    plt.savefig(
        join(pathToSimFolder, f"fitnessDistr_{gate}_pca.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        fitness_bins[:-1],
        fitness_bins,
        weights=fitness_counts,
        color=color(0, 2),
        histtype="step",
        label=r"estimated",
    )
    ax.hist(
        realFitness_bins[:-1],
        realFitness_bins,
        weights=realFitness_counts,
        color=color(1, 2),
        histtype="step",
        label=r"real",
    )

    # ax.set_xlim(0.4,1)
    ax.set_ylim(max(ax.get_ylim()[0], 1e-3), ax.get_ylim()[1])

    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$P(f)$")

    ax.legend(loc="best")

    print(
        f"{gate}: {np.sum(totProbab[hitIndices]):%}    {np.sum(realFitness>minFitness)/samples:%}"
    )

    fitness_counts, fitness_bins = np.histogram(
        estimatedFitness.flatten(),
        weights=totProbab.flatten(),
        bins=fitnessBins,
        density=True,
    )
    realFitness_counts, realFitness_bins = np.histogram(
        realFitness, bins=fitnessBins, density=True
    )

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        fitness_bins[:-1],
        fitness_bins,
        weights=fitness_counts,
        color=color(0, 2),
        histtype="step",
        label=r"estimated",
    )
    ax.hist(
        realFitness_bins[:-1],
        realFitness_bins,
        weights=realFitness_counts,
        color=color(1, 2),
        histtype="step",
        label=r"real",
    )

    # ax.set_xlim(0.4,1)
    # ax.set_ylim(0,1)

    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$P(f)$")

    ax.legend(loc="best")

    plt.savefig(
        join(pathToSimFolder, f"fitnessDistr_{gate}_pca.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.hist(
        fitness_bins[:-1],
        fitness_bins,
        weights=fitness_counts,
        color=color(0, 2),
        histtype="step",
        label=r"estimated",
    )
    ax.hist(
        realFitness_bins[:-1],
        realFitness_bins,
        weights=realFitness_counts,
        color=color(1, 2),
        histtype="step",
        label=r"real",
    )

    # ax.set_xlim(0.4,1)
    ax.set_ylim(max(ax.get_ylim()[0], 1e-3), ax.get_ylim()[1])

    ax.set_xlabel(r"$f$")
    ax.set_ylabel(r"$P(f)$")

    ax.legend(loc="best")

    ax.set_yscale("log")

    plt.savefig(
        join(pathToSimFolder, f"fitnessDistr_{gate}_pca_log.png"),
        bbox_inches="tight",
        dpi=300,
    )
    # plt.show()
    plt.close(fig)

    ### P(f>f_min)

    f = np.linspace(0, 1, 200)

    fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

    ax.plot(
        f,
        [np.sum(totProbab[np.where(estimatedFitness > fi)]) for fi in f],
        color=color(0, 2),
        label=r"estimated",
    )
    ax.plot(
        f,
        [np.sum(realFitness > fi) / samples for fi in f],
        color=color(1, 2),
        label=r"real",
    )

    np.save(
        join(pathToSimFolder, f"{gate}_estimated_integrated_raw_fitness.npy"),
        [np.sum(totProbab[np.where(estimatedFitness > fi)]) for fi in f],
    )
    np.save(
        join(pathToSimFolder, f"{gate}_real_integrated_raw_fitness.npy"),
        [np.sum(realFitness > fi) / samples for fi in f],
    )

    # ax.set_xlim(0.4,1)
    ax.set_ylim(2e-4, None)

    ax.set_xlabel(r"$f_\textrm{min}$")
    ax.set_ylabel(r"$P(f > f_\textrm{min})$")

    ax.grid(which="both")

    ax.legend(loc="best")

    ax.set_yscale("log")

    plt.savefig(
        join(pathToSimFolder, f"P_fmin_{gate}.png"), bbox_inches="tight", dpi=300
    )
    # plt.show()
    plt.close(fig)


############################### raw fitness vs corrected fitness ###############################


rawFitness = np.load(join(pathToSimFolder, "fitness.npy"))
correctedFitness = rawFitness - 2 * np.load(join(pathToSimFolder, "fitnessUncert.npy"))

print(f"corrected fitness prob: {np.sum(correctedFitness>minFitness)/samples:%}")


rawFitness_counts, rawFitness_bins = np.histogram(
    rawFitness, range=(0, 1), bins=fitnessBins, density=True
)
correctedFitness_counts, correctedFitness_bins = np.histogram(
    correctedFitness, range=(0, 1), bins=fitnessBins, density=True
)


fig, ax = plt.subplots(1, 1, figsize=(4.980614173228346, 3.2))

ax.hist(
    rawFitness_bins[:-1],
    rawFitness_bins,
    weights=rawFitness_counts,
    color=color(0, 2),
    histtype="step",
    label=r"raw",
)
ax.hist(
    correctedFitness_bins[:-1],
    correctedFitness_bins,
    weights=correctedFitness_counts,
    color=color(1, 2),
    histtype="step",
    label=r"corrected",
)

ax.set_xlim(0, 1)
# ax.set_ylim(max(ax.get_ylim()[0],1e-3),ax.get_ylim()[1])

ax.set_xlabel(r"$f$")
ax.set_ylabel(r"$P(f)$")

ax.legend(loc="best")

ax.set_yscale("log")

plt.savefig(
    join(pathToSimFolder, f"rawFitnes_vs_correctedFitness.png"),
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


r = np.corrcoef(D_10_01,D_11_00)[0,1]
print(f"D_10_01 vs D_11_00: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_10_01,D_DIFF)[0,1]
print(f"D_10_01 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
r = np.corrcoef(D_11_00,D_DIFF)[0,1]
print(f"D_11_00 vs D_DIFF: r = {r}, T(r) = {r * np.sqrt(n-2)/np.sqrt(1-r**2)}")
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



# D0_counts, D0_bins = np.histogram(D_00_01, bins = bins, density = True)
# D1_counts, D1_bins = np.histogram(D_01_10, bins = bins, density = True)
# Dd_counts, Dd_bins = np.histogram(D_DIFF , bins = bins, density = True)


# D0, D1, Dd = np.meshgrid((D0_bins[1:]+D0_bins[:-1])/2,
#                          (D1_bins[1:]+D1_bins[:-1])/2,
#                          (Dd_bins[1:]+Dd_bins[:-1])/2, indexing = "ij")


# diffs = np.zeros((bins,bins,bins,2,2,2,2))



# diffs[:,:,:,0,0,0,1] = D0
# diffs[:,:,:,0,0,1,0] = D0 + D1
# diffs[:,:,:,0,0,1,1] = 1/3*D1 + 2/3*D0 - Dd
# diffs[:,:,:,0,1,1,0] = D1
# diffs[:,:,:,0,1,1,1] = 1/3*D1 - 1/3*D0 - Dd
# diffs[:,:,:,1,0,1,1] =-2/3*D1 - 1/3*D0 - Dd
# diffs[:,:,:,0,1,0,0] = -diffs[:,:,:,0,0,0,1]
# diffs[:,:,:,1,0,0,0] = -diffs[:,:,:,0,0,1,0]
# diffs[:,:,:,1,1,0,0] = -diffs[:,:,:,0,0,1,1]
# diffs[:,:,:,1,0,0,1] = -diffs[:,:,:,0,1,1,0]
# diffs[:,:,:,1,1,0,1] = -diffs[:,:,:,0,1,1,1]
# diffs[:,:,:,1,1,1,0] = -diffs[:,:,:,1,0,1,1]


# Delta = np.max(np.abs(diffs),axis=(3,4,5,6))

# I_normed = np.zeros((bins,bins,bins,2,2))

# I_normed[:,:,:,0,0] = np.max(diffs[:,:,:,0,0,:,:],axis=(3,4))/Delta
# I_normed[:,:,:,0,1] = np.max(diffs[:,:,:,0,1,:,:],axis=(3,4))/Delta
# I_normed[:,:,:,1,0] = np.max(diffs[:,:,:,1,0,:,:],axis=(3,4))/Delta
# I_normed[:,:,:,1,1] = np.max(diffs[:,:,:,1,1,:,:],axis=(3,4))/Delta


fitness = 1-(I_normed[:,:,:,0,0] + I_normed[:,:,:,0,1] + I_normed[:,:,:,1,0] + 1-I_normed[:,:,:,1,1])/4 #AND




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

ax.legend(loc=2)

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
"""
