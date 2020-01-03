#!/usr/bin/python3

from os.path import join
from re import sub

# class electrodeParameters:
#     def __init__(self,pos,edge,voltage):
#         self.pos    =pos
#         self.edge   =edge
#         self.voltage=voltage


def readParameters(pathToSimFolder):
    parameters={}
    electrodes=[]

    with open(join(pathToSimFolder,"in.txt")) as parameterFile:
        for line in parameterFile:
            line = sub(r' #.*', "", line)
            line = sub(r'#.*', "", line)

            splitted=line.split(" ")
            
            if len(splitted)==2:
                parameters[splitted[0]]=float(splitted[1])
            elif splitted[0]=="electrode":
                electrodes.append([float(splitted[1]),float(splitted[2]),float(splitted[3])])

    return parameters,electrodes



