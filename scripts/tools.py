#!/usr/bin/python3

from os.path import join
from re import sub
import numpy as np
import matplotlib as ma

# ma.use("agg")
import matplotlib.pylab as plt


def readParameters(pathToSimFolder):
    parameters = {}
    electrodes = []

    with open(join(pathToSimFolder, "in.txt")) as parameterFile:
        for line in parameterFile:
            line = sub(r" #.*", "", line)
            line = sub(r"#.*", "", line)

            if line == "\n":
                continue

            splitted = line.split(" ")
            # print(splitted)
            if splitted[0] == "electrode":
                try:
                    geom = parameters["geometry"]
                except KeyError:
                    geom = "rect"
                if geom == "rect":
                    electrodes.append(
                        [float(splitted[1]), float(splitted[2]), float(splitted[3])]
                    )
                elif geom == "circle":
                    electrodes.append([float(splitted[1]), 0, float(splitted[2])])
            else:
                try:
                    parameters[splitted[0]] = float(splitted[1])
                except ValueError:
                    parameters[splitted[0]] = str(splitted[1]).replace("\n", "")
                except IndexError as e:
                    raise ValueError("cant read line: \n\n" + line)

    return parameters, electrodes


def color(i, N):
    if N == 1:
        return ["k"][i]
    if N == 2:
        return ["darkblue", "darkred"][i]
    if N == 3:
        return ["darkgreen", "darkblue", "darkred"][i]
    if N == 4:
        return ["darkgreen", "darkblue", "darkred", "darkorange"][i]
    if N == 5:
        return ["darkgreen", "darkcyan", "darkblue", "darkred", "darkorange"][i]
    if N == 6:
        return ["darkgreen", "darkcyan", "darkblue", "darkred", "darkorange", "gold"][i]
    if N == 7:
        return [
            "darkgreen",
            "darkcyan",
            "darkblue",
            "purple",
            "darkred",
            "darkorange",
            "gold",
        ][i]
    if N == 8:
        return [
            "darkgreen",
            "darkcyan",
            "blue",
            "darkblue",
            "purple",
            "darkred",
            "darkorange",
            "gold",
        ][i]
    if N == 9:
        return [
            "darkgreen",
            "darkcyan",
            "blue",
            "darkblue",
            "purple",
            "darkred",
            "orangered",
            "darkorange",
            "gold",
        ][i]
    if N == 10:
        return [
            "darkgreen",
            "limegreen",
            "darkcyan",
            "blue",
            "darkblue",
            "purple",
            "darkred",
            "orangered",
            "darkorange",
            "gold",
        ][i]

    Nmax = 10
    from matplotlib import colors
    from scipy.interpolate import interp1d

    rgbaData = []
    for j in range(Nmax):
        rgbaData.append(colors.to_rgba(color(j, Nmax)))
    rgbaData = np.array(rgbaData)

    x = np.linspace(0, 1, Nmax)
    fr = interp1d(x, rgbaData[:, 0], kind="cubic")
    fg = interp1d(x, rgbaData[:, 1], kind="cubic")
    fb = interp1d(x, rgbaData[:, 2], kind="cubic")
    RGBA = np.array([fr(i / (N - 1)), fg(i / (N - 1)), fb(i / (N - 1)), 1])
    RGBA[np.where(RGBA > 1)] = 1
    RGBA[np.where(RGBA < 0)] = 0
    return RGBA


def logical(a, b, gate):
    if gate == "AND":
        return np.logical_and(a, b)
    elif gate == "NAND":
        return not np.logical_and(a, b)
    elif gate == "OR":
        return np.logical_or(a, b)
    elif gate == "NOR":
        return not np.logical_or(a, b)
    elif gate == "XOR":
        return np.logical_xor(a, b)
    elif gate == "NXOR":
        return not np.logical_xor(a, b)


from matplotlib.ticker import FuncFormatter


def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(np.around(100 * y))

    # The percent symbol needs escaping in latex
    if ma.rcParams["text.usetex"] is True:
        return s + r"$\%$"
    else:
        return s + "%"


percentFormatter = FuncFormatter(to_percent)
