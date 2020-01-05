#!/usr/bin/python3

from os.path import join
from re import sub
import numpy as np


def readParameters(pathToSimFolder):
    parameters={}
    electrodes=[]

    with open(join(pathToSimFolder,"in.txt")) as parameterFile:
        for line in parameterFile:
            line = sub(r' #.*', "", line)
            line = sub(r'#.*', "", line)

            splitted=line.split(" ")
            
            if len(splitted)==2:
                try:
                    parameters[splitted[0]]=float(splitted[1])
                except ValueError:
                    pass
            elif splitted[0]=="electrode":
                electrodes.append([float(splitted[1]),float(splitted[2]),float(splitted[3])])

    return parameters,electrodes




def color(i,N):
    if N==1: return  ["k"][i]
    if N==2: return  ["darkblue","darkred"][i]
    if N==3: return  ["darkgreen","darkblue","darkred"][i]
    if N==4: return  ["darkgreen","darkblue","darkred","darkorange"][i]
    if N==5: return  ["darkgreen","darkcyan","darkblue","darkred","darkorange"][i]
    if N==6: return  ["darkgreen","darkcyan","darkblue","darkred","darkorange","gold"][i]
    if N==7: return  ["darkgreen","darkcyan","darkblue","purple","darkred","darkorange","gold"][i]
    if N==8: return  ["darkgreen","darkcyan","blue","darkblue","purple","darkred","darkorange","gold"][i]
    if N==9: return  ["darkgreen","darkcyan","blue","darkblue","purple","darkred","orangered","darkorange","gold"][i]
    if N==10:return  ["darkgreen","limegreen","darkcyan","blue","darkblue","purple","darkred","orangered","darkorange","gold"][i]
    
    Nmax=10
    from matplotlib import colors
    from scipy.interpolate import interp1d
    
    rgbaData=[]
    for j in range(Nmax):
        rgbaData.append(colors.to_rgba(color(j,Nmax)))
    rgbaData=np.array(rgbaData)

    x=np.linspace(0,1,Nmax)
    fr=interp1d(x,rgbaData[:,0], kind="cubic")
    fg=interp1d(x,rgbaData[:,1], kind="cubic")
    fb=interp1d(x,rgbaData[:,2], kind="cubic")
    RGBA=np.array([fr(i/(N-1)),fg(i/(N-1)),fb(i/(N-1)),1])
    RGBA[np.where(RGBA>1)]=1
    RGBA[np.where(RGBA<0)]=0
    return RGBA
