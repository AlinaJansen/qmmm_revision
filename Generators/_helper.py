import collections
import datetime

import numpy as np


def _flatten(x):
    """Replace deprecated ``compiler.ast.flatten``"""
    for e in x:
        if not isinstance(e, collections.abc.Iterable) or isinstance(e, str):
            yield e
        else:
            yield from _flatten(e)

def create_dict(label, info):
    out_dict = {}
    for i in range(len(label)):
        out_dict[label[i]] = info[i]
    return out_dict

def logger(log, logstring):
    with open(log, "a") as ofile:
        ofile.write(str(datetime.datetime.now()) + " " + logstring)


def stepper(filename, step):
    """Move to more appropriate module"""
    if step == 0:
        return filename
    elif filename[-4:] == ".g96":
        return str(filename[:-4] + "." + str(step) + ".g96")
    elif filename[-4:] == ".gro":
        return str(filename[:-4] + "." + str(step) + ".gro")
    elif filename[-13:] == ".pointcharges":
        return str(filename[:-13] + "." + str(step) + ".pointcharges")
    else:
        if len(filename) > 7:
            if filename[-7:] == ".fort.7":
                new_filename = str(filename[:-7] + "." + str(step) + ".fort.7")
                return new_filename
        for i in range(0, len(filename)):
            buffer = 0
            if filename[i] == ".":
                buffer = i
        new_filename = str(filename[:buffer] + "." + str(step) + filename[buffer:])
        return new_filename


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

def get_linkatoms_ang(xyzq, qmatomlist, m1list, connlist, prev_scalefacs):
    """Move to more appropriate module"""
    linkatoms = []
    pairlist = []
    for entry in m1list:
        pair = []
        for i in range(0, len(connlist)):
            if int(entry) in np.array(connlist[i]).astype(int):
                if int(entry) == int(connlist[i][0]):
                    for j in range(0, len(connlist[i])):
                        if int(connlist[i][j]) in np.array(qmatomlist).astype(int):
                            pair = [int(connlist[i][j]), int(entry)]
                            pairlist.append(pair)
                            break
                    break
                else:
                    if int(connlist[i][0]) in np.array(qmatomlist).astype(int):
                        pair = [int(connlist[i][0]), int(entry)]
                        pairlist.append(pair)
                        break
    count = 0
    for entry in pairlist:
        linkcoords = []
        q = [
            float(xyzq[entry[0] - 1][0]),
            float(xyzq[entry[0] - 1][1]),
            float(xyzq[entry[0] - 1][2]),
        ]
        m = [
            float(xyzq[entry[1] - 1][0]),
            float(xyzq[entry[1] - 1][1]),
            float(xyzq[entry[1] - 1][2]),
        ]
        qmvec = np.array(m) - np.array(q)
        scalefac = 0.71290813568205  # ratio between B3LYP/6-31G* optimum of C-C in butane vs C-Link relaxed butane

        if len(prev_scalefacs) != 0:
            scalefac = prev_scalefacs[count][3]
        qmvec_norm_scale = np.array(qmvec) * scalefac
        linkcoords = [
            float(qmvec_norm_scale[0]) + float(q[0]),
            float(qmvec_norm_scale[1]) + float(q[1]),
            float(qmvec_norm_scale[2]) + float(q[2]),
            float(scalefac),
        ]
        linkatoms.append(linkcoords)
        count += 1
    return linkatoms

def filter_xyzq(xyzq, atom_list, coordinates=True, charges=True):
    # XX AJ I'm not happy with how this function works (combining lists and arrays) but that's because of different inputs
    python_atom_list = [index - 1 for index in atom_list]
    if isinstance(xyzq, np.ndarray):
        xyzq = xyzq.tolist()
    filtered_xyzq = [xyzq[index] for index in python_atom_list]
    begin = 0 if coordinates else 3
    end = 4 if charges else 3
    filtered_information = np.array(filtered_xyzq)[:, begin:end]
    return filtered_information

def normalized_vector(vec):
    '''
    ------------------------------ \\
    EFFECT: \\
    normalizing a vector \\
    --------------- \\
    INPUT: \\
    vec: np.array
    ------------------------------ \\
    RETURN: \\
    np.array -> normalized vector \\
    --------------- \\
    '''
    return np.array(vec) / np.linalg.norm(np.array(vec))
