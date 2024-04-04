#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA // 
__author__ = 'Alina Jansen'
__date__ = '2024-03-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import sys
import numpy as np

#   Imports From Existing Libraries 

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from Logging.Logger import Logger


#   // TODOS & NOTES //
#   TODO: 
#   - change 'make_xyzq' and 'make_xyzq_io' to be one function
#   NOTE: For now I'm just randomly keeping functions regarding geometry stuff here. Sorting into a class later (AJ)

def readg96(inp):
    '''
    ------------------------------
    EFFECT: \\
    --------------- 
    reads coordinates of all atoms from a g96 file
    ------------------------------
    INPUT: \\
    --------------- 
    inp: string, name of g96 file
    ------------------------------
    RETURN: \\
    --------------- 
    coords: list of coordinates
    ------------------------------
    '''
    coords = []
    with open(inp) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if count == 4:
                break
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                break
    return coords

def readgeo(inp):
    coords = []
    n_a = 0
    with open(inp) as ifile:
        for line in ifile:
            break
        for line in ifile:
            match = re.search(r"^\s*(\d+)", line, flags=re.MULTILINE)
            if match:
                n_a = int(match.group(1))
                break
            else:
                print(".gro is corrupt (no number of atoms found, second line). Exiting.")
                exit(1)
        count = 1
        for line in ifile:
            match = re.search(
                r"^(.{5})(.{5})(.{5})(.{5})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                coords.append(float(match.group(5)) * 10.0)
                coords.append(float(match.group(6)) * 10.0)
                coords.append(float(match.group(7)) * 10.0)
            else:
                print(".gro is corrupt. Exiting.")
                print("Last line:")
                print(line)
                exit(1)
            count += 1
            if count > n_a:
                break
    return coords

def make_xyzq(geo, chargevec):
    xyzq=np.zeros((len(chargevec),4))
    try:
        xyzq[:,0]=geo[0::3]
        xyzq[:,1]=geo[1::3]
        xyzq[:,2]=geo[2::3]
        xyzq[:,3]=chargevec
    except:
        print("Error: Can't make XYZ-charge matrix. Maybe the number of atoms is different in the structure and the topology file ?!")
    return xyzq

def make_xyzq_io(geo, chargevec, outerlist):
    xyzq=np.zeros((len(chargevec),4))
    try:
        xyzq[:,0]=geo[0::3]
        xyzq[:,1]=geo[1::3]
        xyzq[:,2]=geo[2::3]
        outer=np.ones(len(chargevec))
        outer[np.array(outerlist)-1]=0          #the charge of the atoms of the outerlist is set to 0
        xyzq[:,3]=np.array(chargevec)*outer
    except:
        print("Error: Can't make XYZ-charge matrix. Maybe the number of atoms is different in the structure and the topology file ?!")
    return xyzq