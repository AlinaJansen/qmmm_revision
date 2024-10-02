# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

# this will read prepared fields and files for a qmmm job based on gmx.
# will not work alone!

__author__ = "jangoetze"
__date__ = "$06-Feb-2018 12:45:17$"

import math
import os
import re
import subprocess
import copy
import numpy as np
import sqlite3
import sys
import time

from gmx2qmmm._helper import logger, _flatten, stepper
from gmx2qmmm._helper import get_linkatoms_ang, make_xyzq, make_xyzq_io
from gmx2qmmm.operations import expansion_check as rot
from gmx2qmmm.operations import nma_stuff
from gmx2qmmm.operations import nma_3N_6dof as nma
from gmx2qmmm.operations import hes_xyz_g09RevD_01_fchk as write_hess
from gmx2qmmm.operations import bmat                #ness for scan SP
from gmx2qmmm.operations import reindexing # Added by Nicola
from gmx2qmmm.pointcharges import generate_pcf_from_top as make_pcf
from gmx2qmmm.pointcharges import prepare_pcf_for_shift as prep_pcf
from gmx2qmmm.pointcharges import generate_charge_shift as final_pcf

#Ultilites
def get_m2charges(xyzq, m1list, m2list):
    m2charges = []
    count = 0
    for element in m1list:
        m2chargeline = []
        for i in range(0, len(m2list[count])):
            m2chargeline.append(float(xyzq[int(m2list[count][i]) - 1][3]))
        count += 1
        m2charges.append(m2chargeline)
    return m2charges

def read_pcffile(pcffile):
    pcf = []
    with open(pcffile) as ifile:
        for line in ifile:
            match = re.search(r"^QM", line, flags=re.MULTILINE)
            if match:
                pcf.append(["QM"])
                continue
            match = re.search(
                r"^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
            )
            if match:
                pcf.append(
                    [
                        float(match.group(1)),
                        float(match.group(2)),
                        float(match.group(3)),
                        float(match.group(4)),
                    ]
                )
                continue
            match = re.search(r"^$end", line, flags=re.MULTILINE)
            if match:
                break
    return pcf

def read_pcf_self(qmfile):
    pcf_self = 0.0
    with open(qmfile + ".log") as ifile:
        for line in ifile:
            match = re.search(
                r"^\s+Self\s+energy\s+of\s+the\s+charges\s+=\s+([-]*\d+\.\d+)\s+a\.u\.",
                line,
                flags=re.MULTILINE,
            )
            if match:
                pcf_self = float(match.group(1))
                break
    return pcf_self

def remove_inactive(total_force, active):
    new_total_force = []
    for i in range(0, len(total_force)):
        if (i + 1) in np.array(active).astype(int):
            new_total_force.append(total_force[i])
        else:
            new_total_force.append([0.0, 0.0, 0.0])
    return new_total_force

def make_clean_force(total_force):

    clean_force = []
    for element in np.array(total_force):
        forceline = []
        for entry in np.array(element):
            forceline.append(float(entry))
        clean_force.append(forceline)
    return clean_force

def get_full_coords_angstrom(gro):
    fullcoords = []
    with open(gro) as ifile:
        count = 0
        for line in ifile:
            count += 1
            if count == 4:
                break
        count = 0
        for line in ifile:
            count += 1
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if not match:
                break
            full_coord_line = [
                float(match.group(5)) * 10.0,
                float(match.group(6)) * 10.0,
                float(match.group(7)) * 10.0,
            ]
            fullcoords.append(full_coord_line)
    return fullcoords

def get_atoms(qmmmtop, logfile):
    atoms = []
    mass_map = {
        "H": "1.008",
        "He": "4.0026",
        "Li": "6.94",
        "Be": "9.0122",
        "B": "10.81",
        "C": "12.011",
        "N": "14.007",
        "O": "15.999",
        "F": "18.998",
        "Ne": "20.180",
        "Na": "22.990",
        "Mg": "24.305",
        "Al": "26.982",
        "Si": "28.085",
        "P": "30.974",
        "S": "32.06",
        "Cl": "35.45",
        "Ar": "39.948",
        "K": "39.098",
        "Ca": "40.0784",
        "Sc": "44.956",
        "Ti": "47.867",
        "V": "50.942",
        "Cr": "51.996",
        "Mn": "54.938",
        "Fe": "55.8452",
        "Co": "58.933",
        "Ni": "58.693",
        "Cu": "63.5463",
        "Zn": "65.382",
        "Ga": "69.723",
        "Ge": "72.6308",
        "As": "74.922",
        "Se": "78.9718",
        "Br": "79.904",
        "Kr": "83.7982",
        "Rb": "85.468",
        "Sr": "87.62",
        "Y": "88.906",
        "Zr": "91.2242",
        "Nb": "92.906",
        "Mo": "95.95",
        "Tc": "98.906254721",
        "Ru": "101.072",
        "Rh": "102.91",
        "Pd": "106.42",
        "Ag": "107.87",
        "Cd": "112.41",
        "In": "114.82",
        "Sn": "118.71",
        "Sb": "121.76",
        "Te": "127.603",
        "I": "126.90",
        "Xe": "131.29",
        "Cs": "132.91",
        "Ba": "137.33",
        "La": "138.91",
        "Ce": "140.12",
        "Pr": "140.91",
        "Nd": "144.24",
        "Pm": "144.9127493",
        "Sm": "150.362",
        "Eu": "151.96",
        "Gd": "157.253",
        "Tb": "158.93",
        "Dy": "162.50",
        "Ho": "164.93",
        "Er": "167.26",
        "Tm": "168.93",
        "Yb": "173.05",
        "Lu": "174.97",
        "Hf": "178.492",
        "Ta": "180.95",
        "W": "183.84",
        "Re": "186.21",
        "Os": "190.233",
        "Ir": "192.22",
        "Pt": "195.08",
        "Au": "196.97",
        "Hg": "200.59",
        "Tl": "204.38",
        "Pb": "207.2",
        "Bi": "208.98",
        "Po": "208.982430420",
        "At": "209.9871488",
        "Rn": "222.017577725",
        "Fr": "223.019735926",
        "Ra": "226.025409825",
        "Ac": "227.027752126",
        "Th": "232.04",
        "Pa": "231.04",
        "U": "238.03",
        "Np": "237.04817342",
        "Pu": "244.0642045",
        "Am": "243.061381125",
        "Cm": "247.0703545",
        "Bk": "247.0703076",
        "Cf": "251.0795875",
        "Es": "252.082985",
        "Fm": "257.0951067",
        "Md": "258.0984315",
        "No": "259.1010311",
        "Lr": "266.1198356",
        "Rf": "267.1217962",
        "Db": "268.1256757",
        "Sg": "269.1286339",
        "Bh": "270.1333631",
        "Hs": "277.1519058",
        "Mt": "278.1563168",
        "Ds": "281.1645159",
        "Rg": "282.1691272",
        "Cn": "285.177126",
        "Nh": "286.1822172",
        "Fl": "289.190426",
        "Mc": "289.1936389",
        "Lv": "293.204496",
        "Ts": "294.2104674",
        "Og": "295.2162469",
    }
    name_map = {value: key for key, value in mass_map.items()}

    with open(qmmmtop) as ifile:
        for line in ifile:
            match = re.search(r"\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
            if match:
                logger(logfile, "moleculetype section was identified\n")
                break
        for line in ifile:
            match = re.search(r"\[\s+atoms\s*\]", line, flags=re.MULTILINE)
            if match:
                logger(logfile, "atoms section was identified\n")
                break
        for line in ifile:
            match = re.search(r"^\s*\[", line, flags=re.MULTILINE)
            if match:
                break
            match = re.search(
                r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+[\.]*[\d+]*)\s+(\d+[\.]*[\d+]*)",
                line,
                flags=re.MULTILINE,
            )
            if match:
                atomtype = str(match.group(2))
                atommass = float(match.group(8))
                foundname = ""
                # find atom type based on mass
                for key in name_map.items():
                    foundmass = key[0]
                    massdiff = math.sqrt(
                        (float(atommass) - float(foundmass))
                        * (float(atommass) - float(foundmass))
                    )
                    if massdiff < 0.05:
                        foundname = name_map[foundmass]
                        break
                if foundname != "":
                    testmass = mass_map[foundname]
                    massdiff = math.sqrt(
                        (float(atommass) - float(foundmass))
                        * (float(atommass) - float(foundmass))
                    )
                    if massdiff > 0.01:
                        logger(
                            logfile,
                            str(
                                "Found a mass of "
                                + str(atommass)
                                + " for atom type "
                                + str(atomtype)
                                + ' (identified as atom name "'
                                + str(foundname)
                                + '"), which is more than 0.01 different from the expected mass of '
                                + str(testmass)
                                + ". Atom index was "
                                + str(match.group(1))
                                + ". This has no effect except unless the atom was identified wrongly or dynamics are intended. Clean your ffnonbonded.itp to avoid these messages!\n"
                            ),
                        )
                    atoms.append(foundname)
                else:
                    logger(
                        logfile,
                        str(
                            "Atom type "
                            + str(atomtype)
                            + " could not be translated to a regular atom name. Exiting. Last line:\n"
                        ),
                    )
                    logger(logfile, line)
                    exit(1)
    return atoms

def get_nbradius(gro):

    fullcoords = np.array(get_full_coords_nm(gro))
    mindist = fullcoords[0]
    maxdist = fullcoords[1]
    for element in fullcoords:
        for i in range(0, 3):
            if float(mindist[i]) > float(element[i]):
                mindist[i] = element[i]
            if float(maxdist[i]) < float(element[i]):
                maxdist[i] = element[i]
    maxcoords = np.array(maxdist) - np.array(mindist)
    return np.linalg.norm(maxcoords)

def get_full_coords_nm(gro):  # read g96
    fullcoords = []
    with open(gro) as ifile:
        count = 0
        for line in ifile:
            if count == 4:
                break
            count += 1
        for line in ifile:
            match = re.search(
                r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                line,
                flags=re.MULTILINE,
            )
            if not match:
                break
            full_coord_line = [
                float(match.group(5)),
                float(match.group(6)),
                float(match.group(7)),
            ]
            fullcoords.append(full_coord_line)
    return fullcoords

def get_approx_hessian(xyz, old_xyz, grad, old_grad, old_hess, logfile):
    s = xyz - old_xyz
    if s.shape[1] != 1:
        s = s.reshape(3 * len(s), 1)  # reshape s so it can be used in dot products
    g = grad - old_grad
    if g.shape[1] != 1:
        g = g.reshape(3 * len(g), 1)  # reshape g so it can be used in dot products
    # Broyden-Fletcher-Goldfarb-Shanno
    # Define a few values and matrices first for convenience
    mat = old_hess.dot(s)
    factor1 = g.T.dot(s)
    factor2 = s.T.dot(mat)
    idmat = np.eye(len(g))

    # update formula
    new_hess = old_hess + g.dot(g.T) / factor1 - mat.dot(mat.T) / factor2
    # moreover, we calculate eigenvalues as they are indicative of the curvature of the current PES
    eigvals, eigvecs = np.linalg.eig(new_hess)

    # Check BFGS condition
    if factor1 > 0:
        logger(logfile, "BFGS condition fulfilled.\n")
        WARN = False
    else:
        logger(logfile, "BFGS condition not fulfilled! We keep the old Hessian.\n")
        WARN = True
        new_hess = old_hess

    return new_hess, eigvals.min(), WARN

def calculate_nuclear_self_interaction(inp):
    import numpy as np
    xyzq = []
    e = 0
    with open(inp) as ifile:
        lines = ifile.readlines()
        for line in lines:
            data = line.strip().split()
            if len(data) == 5:
                if data[0] == 'Q':
                    try:
                        x, y, z = float(data[2]), float(data[3]), float(data[4])
                        q = float(data[1])
                        xyzq.append([x, y, z, q])
                    except:
                        continue
                else:
                    continue
            else:
                continue

        xyzarr = np.array(xyzq).astype(float)
        for i in range(len(xyzarr)):
            for j in range(i+1, len(xyzarr), 1):
                dist = bmat.length(xyzarr[i][:3], xyzarr[j][:3])
                qprod = xyzarr[i][3] * xyzarr[j][3]
                col = qprod/dist if dist != 0 else 0
                e += col

    return e

def get_working_directory():
    return os.getcwd()

def extract_qmsubsystem_from_qmlist(qmlist, a, b):
    endqmlist = []
    foundbegin = False
    counter = 0
    
    for i in range(len(qmlist)):
        if qmlist[i] == a:
            foundbegin = True
            counter = i
            break
    
    if foundbegin:
        for i in range(counter, len(qmlist)):
            endqmlist.append(qmlist[i])
            if qmlist[i] == b:
                break

    return endqmlist

def pipein_define_TM(jobname, inpbasisset, inpcharge, inpgeneralmethod, inpfunctional):
    '''commands = f"""
    {jobname}
    aa str
    ired
    *
    b all {inpbasisset}
    *
    eht

    {inpcharge}

    {inpgeneralmethod}
    on
    func {inpfunctional}
    *
    ri
    on
    *
    *
    """'''
    commands = f"""
    {jobname}
    aa str
    ired
    *
    b all {inpbasisset}
    *
    eht

    {inpcharge}

    {inpgeneralmethod}
    on
    func {inpfunctional}
    *
    ri
    on
    *
    *
    """

    # Start the `define` process
    process = subprocess.Popen(
        ['define'],  # Start `define` program
        stdin=subprocess.PIPE,  # Enable writing input to the process
        stdout=subprocess.PIPE,  # Capture the output
        stderr=subprocess.PIPE,  # Capture any errors
        text=True,  # Treat input/output as text instead of bytes
    )

    # Send the commands to the `define` process
    stdout, stderr = process.communicate(input=commands)

    print(stdout)
    #if stderr:
    #    print(f"Error: {stderr}")

def make_TM_inp(qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    #qmmmtop = qmmmInputs.top                       SP activate next line
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag
    connlist = qmmmInputs.connlist
    xyzq = qmmmInputs.xyzq
    m1list = qmmmInputs.m1list

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra
    generalmethod = qmmmInputs.qmparams.generalmethod
    maxcyc = qmmmInputs.qmparams.maxcyc

    insert = ""
    oldinsert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step) ))
        if int(curr_step) > 1:
            oldinsert = str("." + str(int(curr_step) - 1))
    tmfile = "control"
    mosfile = str(jobname + insert + ".mos")
    oldmosfile = str(jobname + oldinsert + ".mos")

    fullcoords = get_full_coords_angstrom(gro)
    atoms = get_atoms(qmmmtop, logfile)
    original_stdout = sys.stdout # Save a reference to the original standard output
    with open('top_info.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print('QMMMTOPINFO')
        for i in atoms:
            print(str(i))
        sys.stdout = original_stdout # Reset the standard output to its original value
    with open('gro_info.txt', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print('GROINFO')
        for i in fullcoords:
            print(str(i))
        sys.stdout = original_stdout # Reset the standard output to its original value
    with open("str", 'w') as struc:
        struc.write("$coord\n")
        count = 0
        for element in fullcoords:
            if int(count + 1) in np.array(qmatomlist).astype(int):
                #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                struc.write(
                    "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                        str(atoms[count]),
                        float(element[0]),
                        float(element[1]),
                        float(element[2]),
                    )
                )
            count += 1
        for element in linkatoms:
            # links are already in angstrom
            struc.write(
                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                    str("H"), float(element[0]), float(element[1]), float(element[2])
                )
            )
        struc.write("$end\n")
    time.sleep(10) 
    pipein_define_TM(str(jobname + insert), str(basis), str(charge), str(generalmethod), str(method)) 
    time.sleep(15)  
    if int(curr_step) != 0 or nmaflag == 1:
        subprocess.call("rm mos", shell=True)
    with open(tmfile, 'r') as ifile:
        lines = ifile.readlines()
    with open(tmfile, 'w') as cfile:
        for line in lines:
            if '$scfmo' in line and (int(curr_step) != 0 or nmaflag == 1):
                cfile.write('$scfmo   file=' + oldmosfile + '\n')
                continue
            elif '$scfiterlimit' in line:
                cfile.write("$scfiterlimit       " + str(maxcyc) + "\n")
                continue
            elif '$maxcore' in line:
                cfile.write("$maxcor    " + str(memory) + " MiB  per_core\n")
                continue
            elif '$drvopt' in line:
                cfile.write(line)
                cfile.write("   point charges\n")
                continue
            elif '$grad' in line:
                cfile.write(line)
                cfile.write("$point_charges nocheck list\n")
                with open(pcffile) as ifile:
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            cfile.write(
                                " {:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                    float(match.group(1)),
                                    float(match.group(2)),
                                    float(match.group(3)),
                                    float(match.group(4)),
                                )
                            )
                cfile.write("$mvd\n")
                cfile.write("$moments\n")
                cfile.write("$pop\n")
                cfile.write("$esp_fit kollman\n")
                continue
            elif '$scfconv' in line:
                cfile.write("$scfconv   8\n")
                continue
            elif '$ricore' in line:
                cfile.write("$ricore      " + str(memory) + "\n")
                cfile.write(str(extra) + "\n")
                continue
            else:
                cfile.write(line)
    return str(jobname + insert + ".inp")
            


def make_serenity_inp(qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    #qmmmtop = qmmmInputs.top                       SP activate next line
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag
    connlist = qmmmInputs.connlist
    xyzq = qmmmInputs.xyzq
    m1list = qmmmInputs.m1list

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra
    generalmethod = qmmmInputs.qmparams.generalmethod
    maxcyc = qmmmInputs.qmparams.maxcyc
    fde = qmmmInputs.qmparams.fde
    tddft = qmmmInputs.qmparams.tddft
    tddftstates = qmmmInputs.qmparams.tddftstates

    insert = ""
    oldinsert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step) ))
        if int(curr_step) > 1:
            oldinsert = str("." + str(int(curr_step) - 1))
    serenityfile = str(jobname + insert + ".inp")
    #gbwfile = str(jobname + insert + ".gbw")
    #oldgbwfile = str(jobname + oldinsert + ".gbw")

    #if nmaflag == 1:
    #    oldgbwfile = str(jobname + ".gbw")

    with open(serenityfile, "w") as ofile:
        fullcoords = get_full_coords_angstrom(gro)
        atoms = get_atoms(qmmmtop, logfile)
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('top_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('QMMMTOPINFO')
            for i in atoms:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        with open('gro_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('GROINFO')
            for i in fullcoords:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        if not fde:
            ofile.write("+system" + "\n")
            ofile.write("   name " + str(jobname + insert) + "\n")
            #if int(curr_step) != 0 or nmaflag == 1:
            #    wd = get_working_directory()
            #    ofile.write("   load " + wd + "/" + "\n")
            ofile.write("   geometry ./" + str(jobname + insert)+".xyz" + "\n")
        
            ofile.write("   method " + str(generalmethod) + "\n")
            if generalmethod == "dft":
                ofile.write("   +dft" + "\n")
                ofile.write("       functional " + str(method) + "\n")
                ofile.write("   -dft" + "\n")
            ofile.write("   charge " + str(int(charge)) + "\n")
            ofile.write("   spin " + str(float(int(multi - 1)/2)) + "\n")
            ofile.write("   +basis" + "\n")
            ofile.write("       label " + str(basis) + "\n")
            ofile.write("   -basis" + "\n")
            ofile.write("   +grid" + "\n")
            ofile.write("       accuracy 6" + "\n")
            ofile.write("       smallGridAccuracy 3" + "\n")
            ofile.write("   -grid" + "\n")
        
            ofile.write("   +extcharges" + "\n")
            ofile.write("       externalChargesFile pointcharges.pc" + "\n")
            ofile.write("   -extcharges" + "\n")
        
            ofile.write("   +scf" + "\n")
            #ofile.write("       maxCycles " + str(maxcyc) + "\n")
            ofile.write("       initialguess EHT" + "\n")
            #ofile.write("       energyThreshold 1e-6" + "\n")
            ofile.write("   -scf" + "\n")
            ofile.write("-system" + "\n")

            #ofile.write("+task scf" + "\n")
            #ofile.write("   system " + str(jobname + insert) + "\n")
            #ofile.write("-task" + "\n")
            if tddft == 1:
                ofile.write("+task lrscf" + "\n")
                ofile.write("   act " + str(jobname + insert) + "\n")
                ofile.write("   nEigen " + str(tddftstates) + "\n")
                ofile.write("-task" + "\n")

            ofile.write("+task GradientTask" + "\n")
            ofile.write("   system " + str(jobname + insert) + "\n")
            ofile.write("-task" + "\n")

            ofile.write("+task pop" + "\n")
            ofile.write("   act " + str(jobname + insert) + "\n")
            ofile.write("   algorithm Hirshfeld" + "\n")
            ofile.write("-task" + "\n")


            #ofile.write("!" + str(method))
            #if str(basis) != "NONE":
            #    ofile.write(" " + str(basis))
            #ofile.write(" ENGRAD PrintMOs Printbasis NoUseSym SCFCONV6" + "\n")
            #if str(extra) != "NONE":
            #    ofile.write(str(extra) + "\n")
            #if int(curr_step) == 0:
            #    ofile.write("!HUECKEL" + "\n")
            #if int(curr_step) != 0 or nmaflag == 1: #Ask Denis
            #    ofile.write("!MORead" + "\n")
            #    ofile.write("%moinp " + "\"" + oldgbwfile + "\"" +"\n")
        
            #ofile.write("%maxcore" + " " + str(memory) + "\n")
            #ofile.write("%pal" + "\n" + "nprocs" + " " + str(cores) + "\n" + "end" + "\n")
        

        
            #ofile.write("%pointcharges \"pointcharges.pc\" " + "\n\n")
            #ofile.write(
            #    "\n#QM part of step "+str(curr_step)+"\n\n"
            #    + "*"
            #    + " "
            #    + "xyz"
            #    + " "
            #    + str(int(charge))
            #    + " "
            #    + str(int(multi))
            #    + "\n"
            #)

            with open(str(jobname + insert)+".xyz", 'w') as xyzfile:
            
                count1 = 0
                count2 = 0
                for element in fullcoords:
                    if int(count1 + 1) in np.array(qmatomlist).astype(int):
                        count2 += 1
                    count1 +=1
            
                for element in linkatoms:
                    count2 += 1
            
                xyzfile.write(str(count2) + "\n")


                xyzfile.write("Hi there, i am using gmx2qmmm/serenity" + "\n")
                count = 0
                for element in fullcoords:
                    if int(count + 1) in np.array(qmatomlist).astype(int):
                        #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                        xyzfile.write(
                            "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                str(atoms[count]),
                                float(element[0]),
                                float(element[1]),
                                float(element[2]),
                            )
                        )
                    count += 1
                for element in linkatoms:
                    # links are already in angstrom
                    xyzfile.write(
                        "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            str("H"), float(element[0]), float(element[1]), float(element[2])
                        )
                    )
                #ofile.write("\n")
                #ofile.write("*" + "\n")


            with open("pointcharges.pc", 'w') as pfile:

                with open(pcffile) as ifile:
                    c = 0
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            c += 1
            
                with open(pcffile) as ifile:
                    pfile.write(str(int(c)) + "\n")
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:       
                            pfile.write(
                                "{:>12.6f} {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                                    float(match.group(4)),
                                    float(match.group(1)),
                                    float(match.group(2)),
                                    float(match.group(3)),
                                )
                            )
            #ofile.write("\n")
        else:
            #gmxplus = qmmmInputs.mmparams.gmxplus
            inout = qmmmInputs.inout
            #if gmxplus and inout:
            if inout:
                inner = qmmmInputs.innerlistorig
                outer = qmmmInputs.outerlistorig
                qm = qmmmInputs.qmlistorig
            
                ref = qm + inner + outer
                ref = list(set(ref))
                memory_dict = reindexing.reindexing_memory(ref, outer)
            '''ofile.write("+system" + "\n")
            ofile.write("   name " + str(jobname + insert) + "\n")
            #if int(curr_step) != 0 or nmaflag == 1:
            #    wd = get_working_directory()
            #    ofile.write("   load " + wd + "/" + "\n")
            ofile.write("   geometry ./" + str(jobname + insert)+".xyz" + "\n")
        
            ofile.write("   method " + str(generalmethod) + "\n")
            if generalmethod == "dft":
                ofile.write("   +dft" + "\n")
                ofile.write("       functional " + str(method) + "\n")
                ofile.write("   -dft" + "\n")
            ofile.write("   charge " + str(int(charge)) + "\n")
            ofile.write("   spin " + str(float(int(multi - 1)/2)) + "\n")
            ofile.write("   +basis" + "\n")
            ofile.write("       label " + str(basis) + "\n")
            ofile.write("   -basis" + "\n")
            ofile.write("   +grid" + "\n")
            ofile.write("       accuracy 6" + "\n")
            ofile.write("       smallGridAccuracy 3" + "\n")
            ofile.write("   -grid" + "\n")
        
            ofile.write("   +extcharges" + "\n")
            ofile.write("       externalChargesFile pointcharges.pc" + "\n")
            ofile.write("   -extcharges" + "\n")
        
            ofile.write("   +scf" + "\n")
            #ofile.write("       maxCycles " + str(maxcyc) + "\n")
            ofile.write("       initialguess EHT" + "\n")
            #ofile.write("       energyThreshold 1e-6" + "\n")
            ofile.write("   -scf" + "\n")
            ofile.write("-system" + "\n")

            #ofile.write("+task scf" + "\n")
            #ofile.write("   system " + str(jobname + insert) + "\n")
            #ofile.write("-task" + "\n")

            #ofile.write("+task GradientTask" + "\n")
            #ofile.write("   system " + str(jobname + insert) + "\n")
            #ofile.write("-task" + "\n")

            #ofile.write("+task pop" + "\n")
            #ofile.write("   act " + str(jobname + insert) + "\n")
            #ofile.write("   algorithm Hirshfeld" + "\n")
            #ofile.write("-task" + "\n")


            #ofile.write("!" + str(method))
            #if str(basis) != "NONE":
            #    ofile.write(" " + str(basis))
            #ofile.write(" ENGRAD PrintMOs Printbasis NoUseSym SCFCONV6" + "\n")
            #if str(extra) != "NONE":
            #    ofile.write(str(extra) + "\n")
            #if int(curr_step) == 0:
            #    ofile.write("!HUECKEL" + "\n")
            #if int(curr_step) != 0 or nmaflag == 1: #Ask Denis
            #    ofile.write("!MORead" + "\n")
            #    ofile.write("%moinp " + "\"" + oldgbwfile + "\"" +"\n")
        
            #ofile.write("%maxcore" + " " + str(memory) + "\n")
            #ofile.write("%pal" + "\n" + "nprocs" + " " + str(cores) + "\n" + "end" + "\n")
        

        
            #ofile.write("%pointcharges \"pointcharges.pc\" " + "\n\n")
            #ofile.write(
            #    "\n#QM part of step "+str(curr_step)+"\n\n"
            #    + "*"
            #    + " "
            #    + "xyz"
            #    + " "
            #    + str(int(charge))
            #    + " "
            #    + str(int(multi))
            #    + "\n"
            #)

            with open(str(jobname + insert)+".xyz", 'w') as xyzfile:
            
                count1 = 0
                count2 = 0
                for element in fullcoords:
                    if int(count1 + 1) in np.array(qmatomlist).astype(int):
                        count2 += 1
                    count1 +=1
            
                for element in linkatoms:
                    count2 += 1
            
                xyzfile.write(str(count2) + "\n")


                xyzfile.write("Hi there, i am using gmx2qmmm/serenity" + "\n")
                count = 0
                for element in fullcoords:
                    if int(count + 1) in np.array(qmatomlist).astype(int):
                        #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                        xyzfile.write(
                            "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                str(atoms[count]),
                                float(element[0]),
                                float(element[1]),
                                float(element[2]),
                            )
                        )
                    count += 1
                for element in linkatoms:
                    # links are already in angstrom
                    xyzfile.write(
                        "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            str("H"), float(element[0]), float(element[1]), float(element[2])
                        )
                    )'''
                #ofile.write("\n")
                #ofile.write("*" + "\n")


            with open("pointcharges.pc", 'w') as pfile:

                with open(pcffile) as ifile:
                    c = 0
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            c += 1
            
                with open(pcffile) as ifile:
                    pfile.write(str(int(c)) + "\n")
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:       
                            pfile.write(
                                "{:>12.6f} {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                                    float(match.group(4)),
                                    float(match.group(1)),
                                    float(match.group(2)),
                                    float(match.group(3)),
                                )
                            )
            #ofile.write("\n")
            naddxcfunc = qmmmInputs.qmparams.naddxcfunc
            naddkinfunc = qmmmInputs.qmparams.naddkinfunc
            tddftnaddxcfunc = qmmmInputs.qmparams.tddftnaddxcfunc
            tddftnaddkinfunc = qmmmInputs.qmparams.tddftnaddkinfunc
            numactsys = qmmmInputs.qmparams.numactsys
            numenvsys = qmmmInputs.qmparams.numenvsys
            namesact = []
            namesenv = []  
            for i in range(numactsys):
                endactsys = getattr(qmmmInputs.qmparams, f"endactsys{i}")
                generalmethod = getattr(qmmmInputs.qmparams, f"generalmethod{i}")
                maxcyc = getattr(qmmmInputs.qmparams, f"maxcyc{i}")
                method = getattr(qmmmInputs.qmparams, f"method{i}")
                basis = getattr(qmmmInputs.qmparams, f"basis{i}")
                charge = getattr(qmmmInputs.qmparams, f"charge{i}")
                multi = getattr(qmmmInputs.qmparams, f"multiplicity{i}")
                #if gmxplus and inout:
                if inout:
                    qmsublist = extract_qmsubsystem_from_qmlist(qmatomlist, memory_dict[endactsys[0]], memory_dict[endactsys[1]])
                else:
                    qmsublist = extract_qmsubsystem_from_qmlist(qmatomlist, endactsys[0], endactsys[1])
                sublinkatoms = get_linkatoms_ang(xyzq, qmsublist, m1list, connlist, [])

                with open(str(jobname + insert)+ str(i)+".xyz", 'w') as xyzfile:
            
                    count1 = 0
                    count2 = 0
                    for element in fullcoords:
                        if int(count1 + 1) in np.array(qmsublist).astype(int):
                            count2 += 1
                        count1 +=1
            
                    for element in sublinkatoms:
                        count2 += 1
            
                    xyzfile.write(str(count2) + "\n")


                    xyzfile.write("Hi there, i am using gmx2qmmm/serenity" + "\n")
                    count = 0
                    for element in fullcoords:
                        if int(count + 1) in np.array(qmsublist).astype(int):
                            xyzfile.write(
                                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                    str(atoms[count]),
                                    float(element[0]),
                                    float(element[1]),
                                    float(element[2]),
                                )
                            )
                        count += 1
                    for element in sublinkatoms:
                        # links are already in angstrom
                        xyzfile.write(
                            "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                str("H"), float(element[0]), float(element[1]), float(element[2])
                            )
                        )
             
                ofile.write("+system" + "\n")
                ofile.write("   name " + str(jobname + insert) + str(i) + "\n")
                namesact.append(str(jobname + insert + str(i)))
                #if int(curr_step) != 0 or nmaflag == 1:
                #    wd = get_working_directory()
                #    ofile.write("   load " + wd + "/" + "\n")
                ofile.write("   geometry ./" + str(jobname + insert)+ str(i)+".xyz" + "\n")
        
                ofile.write("   method " + str(generalmethod) + "\n")
                if generalmethod == "dft":
                    ofile.write("   +dft" + "\n")
                    ofile.write("       functional " + str(method) + "\n")
                    ofile.write("   -dft" + "\n")
                ofile.write("   charge " + str(int(charge)) + "\n")
                ofile.write("   spin " + str(float(int(multi - 1)/2)) + "\n")
                ofile.write("   +basis" + "\n")
                ofile.write("       label " + str(basis) + "\n")
                ofile.write("   -basis" + "\n")
                ofile.write("   +grid" + "\n")
                ofile.write("       accuracy 6" + "\n")
                ofile.write("       smallGridAccuracy 3" + "\n")
                ofile.write("   -grid" + "\n")
        
                ofile.write("   +extcharges" + "\n")
                ofile.write("       externalChargesFile pointcharges.pc" + "\n")
                ofile.write("   -extcharges" + "\n")
        
                ofile.write("   +scf" + "\n")
                #ofile.write("       maxCycles " + str(maxcyc) + "\n")
                ofile.write("       initialguess EHT" + "\n")
                #ofile.write("       energyThreshold 1e-6" + "\n")
                ofile.write("   -scf" + "\n")
                ofile.write("-system" + "\n")
                ofile.write("\n")

            for i in range(numenvsys):
                endenvsys = getattr(qmmmInputs.qmparams, f"endenvsys{i}")
                generalmethod = getattr(qmmmInputs.qmparams, f"egeneralmethod{i}")
                maxcyc = getattr(qmmmInputs.qmparams, f"emaxcyc{i}")
                method = getattr(qmmmInputs.qmparams, f"emethod{i}")
                basis = getattr(qmmmInputs.qmparams, f"ebasis{i}")
                charge = getattr(qmmmInputs.qmparams, f"echarge{i}")
                multi = getattr(qmmmInputs.qmparams, f"emultiplicity{i}")
                #if gmxplus and inout:
                if inout:
                    qmsublist = extract_qmsubsystem_from_qmlist(qmatomlist, memory_dict[endenvsys[0]], memory_dict[endenvsys[1]])
                else:
                    qmsublist = extract_qmsubsystem_from_qmlist(qmatomlist, endenvsys[0], endenvsys[1])
                sublinkatoms = get_linkatoms_ang(xyzq, qmsublist, m1list, connlist, [])

                with open("e" + str(jobname + insert)+ str(i)+".xyz", 'w') as xyzfile:
            
                    count1 = 0
                    count2 = 0
                    for element in fullcoords:
                        if int(count1 + 1) in np.array(qmsublist).astype(int):
                            count2 += 1
                        count1 +=1
            
                    for element in sublinkatoms:
                        count2 += 1
            
                    xyzfile.write(str(count2) + "\n")


                    xyzfile.write("Hi there, i am using gmx2qmmm/serenity" + "\n")
                    count = 0
                    for element in fullcoords:
                        if int(count + 1) in np.array(qmsublist).astype(int):
                            xyzfile.write(
                                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                    str(atoms[count]),
                                    float(element[0]),
                                    float(element[1]),
                                    float(element[2]),
                                )
                            )
                        count += 1
                    for element in sublinkatoms:
                        # links are already in angstrom
                        xyzfile.write(
                            "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                str("H"), float(element[0]), float(element[1]), float(element[2])
                            )
                        )

                ofile.write("+system" + "\n")
                ofile.write("   name " + "e" + str(jobname + insert) + str(i) + "\n")
                namesenv.append(str("e" + jobname + insert + str(i)))
                #if int(curr_step) != 0 or nmaflag == 1:
                #    wd = get_working_directory()
                #    ofile.write("   load " + wd + "/" + "\n")
                ofile.write("   geometry ./" + "e" + str(jobname + insert)+ str(i)+".xyz" + "\n")
        
                ofile.write("   method " + str(generalmethod) + "\n")
                if generalmethod == "dft":
                    ofile.write("   +dft" + "\n")
                    ofile.write("       functional " + str(method) + "\n")
                    ofile.write("   -dft" + "\n")
                ofile.write("   charge " + str(int(charge)) + "\n")
                ofile.write("   spin " + str(float(int(multi - 1)/2)) + "\n")
                ofile.write("   +basis" + "\n")
                ofile.write("       label " + str(basis) + "\n")
                ofile.write("   -basis" + "\n")
                ofile.write("   +grid" + "\n")
                ofile.write("       accuracy 6" + "\n")
                ofile.write("       smallGridAccuracy 3" + "\n")
                ofile.write("   -grid" + "\n")
        
                ofile.write("   +extcharges" + "\n")
                ofile.write("       externalChargesFile pointcharges.pc" + "\n")
                ofile.write("   -extcharges" + "\n")
        
                ofile.write("   +scf" + "\n")
                #ofile.write("       maxCycles " + str(maxcyc) + "\n")
                ofile.write("       initialguess EHT" + "\n")
                #ofile.write("       energyThreshold 1e-6" + "\n")
                ofile.write("   -scf" + "\n")
                ofile.write("-system" + "\n")
                ofile.write("\n")
            if tddft == 1:
                ofile.write("+task FaT" + "\n")
                for i in namesact:
                    ofile.write("   act " + i + "\n")
                for i in namesenv:
                    ofile.write("   env " + i + "\n")
                ofile.write("   +emb" + "\n")
                ofile.write("       naddXCFunc " + str(tddftnaddxcfunc) + "\n")
                ofile.write("       naddKinFunc " + str(tddftnaddkinfunc) + "\n")
                ofile.write("   -emb" + "\n")
                ofile.write("-task" + "\n")

                for i in namesact:
                    ofile.write("+task lrscf" + "\n")
                    ofile.write("   act " + i + "\n")
                    container = namesact.copy()
                    container.remove(i)
                    for i in container:
                        ofile.write("   env " + i + "\n")
                    for i in namesenv:
                        ofile.write("   env " + i + "\n")
                    ofile.write("   nEigen " + str(tddftstates) + "\n")
                    ofile.write("   +emb" + "\n")
                    ofile.write("       naddXCFunc " + str(tddftnaddxcfunc) + "\n")
                    ofile.write("       naddKinFunc " + str(tddftnaddkinfunc) + "\n")
                    ofile.write("   -emb" + "\n")
                    ofile.write("-task" + "\n")

            ofile.write("+task GradientTask" + "\n")
            for i in namesact:
                ofile.write("   act " + i + "\n")
            for i in namesenv:
                ofile.write("   env " + i + "\n")
            ofile.write("   +emb" + "\n")
            ofile.write("       naddXCFunc " + str(naddxcfunc) + "\n")
            ofile.write("       naddKinFunc " + str(naddkinfunc) + "\n")
            ofile.write("   -emb" + "\n")
            ofile.write("-task" + "\n")
            for i in namesact:
                ofile.write("+task pop" + "\n")
            #ofile.write("   act " + str(jobname + insert) + "\n")

                ofile.write("   act " + i + "\n")

                ofile.write("   algorithm Hirshfeld" + "\n")
                ofile.write("-task" + "\n")
            for i in namesenv:
                ofile.write("+task pop" + "\n")
                ofile.write("   act " + i + "\n")
                ofile.write("   algorithm Hirshfeld" + "\n")
                ofile.write("-task" + "\n")
            
            '''with open("pointcharges.pc", 'w') as pfile:

                with open(pcffile) as ifile:
                    c = 0
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            c += 1
            
                with open(pcffile) as ifile:
                    pfile.write(str(int(c)) + "\n")
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:       
                            pfile.write(
                                "{:>12.6f} {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                                    float(match.group(4)),
                                    float(match.group(1)),
                                    float(match.group(2)),
                                    float(match.group(3)),
                                )
                            )'''

                
        

    return serenityfile
#qm & mm program ultis
'''
def make_serenity_inp2(qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    #qmmmtop = qmmmInputs.top                       SP activate next line
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra
    generalmethod = qmmmInputs.qmparams.generalmethod
    maxcyc = qmmmInputs.qmparams.maxcyc

    insert = ""
    oldinsert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step) ))
        if int(curr_step) > 1:
            oldinsert = str("." + str(int(curr_step) - 1))
    serenityfile = str(jobname + insert + ".inp")
    #gbwfile = str(jobname + insert + ".gbw")
    #oldgbwfile = str(jobname + oldinsert + ".gbw")

    #if nmaflag == 1:
    #    oldgbwfile = str(jobname + ".gbw")

    with open(serenityfile, "w") as ofile:
        fullcoords = get_full_coords_angstrom(gro)
        atoms = get_atoms(qmmmtop, logfile)
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('top_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('QMMMTOPINFO')
            for i in atoms:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        with open('gro_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('GROINFO')
            for i in fullcoords:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        
        ofile.write("+system" + "\n")
        ofile.write("   name " + str(jobname + insert) + "\n")
        #if int(curr_step) != 0 or nmaflag == 1:
        #    wd = get_working_directory()
        #    ofile.write("   load " + wd + "/" + "\n")
        ofile.write("   geometry ./" + str(jobname + insert)+".xyz" + "\n")
        
        ofile.write("   method " + str(generalmethod) + "\n")
        if generalmethod == "dft":
            ofile.write("   +dft" + "\n")
            ofile.write("       functional " + str(method) + "\n")
            ofile.write("   -dft" + "\n")
        ofile.write("   charge " + str(int(charge)) + "\n")
        ofile.write("   spin " + str(float(int(multi - 1)/2)) + "\n")
        ofile.write("   +basis" + "\n")
        ofile.write("       label " + str(basis) + "\n")
        ofile.write("   -basis" + "\n")
        ofile.write("   +grid" + "\n")
        ofile.write("       accuracy 6" + "\n")
        ofile.write("       smallGridAccuracy 3" + "\n")
        ofile.write("   -grid" + "\n")
        
        ofile.write("   +extcharges" + "\n")
        ofile.write("       externalChargesFile pointcharges.pc" + "\n")
        ofile.write("   -extcharges" + "\n")
        
        ofile.write("   +scf" + "\n")
        ofile.write("       maxCycles " + str(maxcyc) + "\n")
        if int(curr_step) == 0:
            ofile.write("       initialguess EHT" + "\n")
        ofile.write("       energyThreshold 1e-6" + "\n")
        ofile.write("   -scf" + "\n")
        ofile.write("-system" + "\n")

        #ofile.write("+task scf" + "\n")
        #ofile.write("   system " + str(jobname + insert) + "\n")
        #ofile.write("-task" + "\n")

        ofile.write("+task GradientTask" + "\n")
        ofile.write("   system " + str(jobname + insert) + "\n")
        ofile.write("-task" + "\n")

        ofile.write("+task pop" + "\n")
        ofile.write("   act " + str(jobname + insert) + "\n")
        ofile.write("   algorithm Hirshfeld" + "\n")
        ofile.write("-task" + "\n")


        #ofile.write("!" + str(method))
        #if str(basis) != "NONE":
        #    ofile.write(" " + str(basis))
        #ofile.write(" ENGRAD PrintMOs Printbasis NoUseSym SCFCONV6" + "\n")
        #if str(extra) != "NONE":
        #    ofile.write(str(extra) + "\n")
        #if int(curr_step) == 0:
        #    ofile.write("!HUECKEL" + "\n")
        #if int(curr_step) != 0 or nmaflag == 1: #Ask Denis
        #    ofile.write("!MORead" + "\n")
        #    ofile.write("%moinp " + "\"" + oldgbwfile + "\"" +"\n")
        
        #ofile.write("%maxcore" + " " + str(memory) + "\n")
        #ofile.write("%pal" + "\n" + "nprocs" + " " + str(cores) + "\n" + "end" + "\n")
        

        
        #ofile.write("%pointcharges \"pointcharges.pc\" " + "\n\n")
        #ofile.write(
        #    "\n#QM part of step "+str(curr_step)+"\n\n"
        #    + "*"
        #    + " "
        #    + "xyz"
        #    + " "
        #    + str(int(charge))
        #    + " "
        #    + str(int(multi))
        #    + "\n"
        #)

        with open(str(jobname + insert)+".xyz", 'w') as xyzfile:
            
            count1 = 0
            count2 = 0
            for element in fullcoords:
                if int(count1 + 1) in np.array(qmatomlist).astype(int):
                    count2 += 1
                count1 +=1
            
            for element in linkatoms:
                count2 += 1
            
            xyzfile.write(str(count2) + "\n")


            xyzfile.write("Hi there, i am using gmx2qmmm/serenity" + "\n")
            count = 0
            for element in fullcoords:
                if int(count + 1) in np.array(qmatomlist).astype(int):
                    #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                    xyzfile.write(
                        "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            str(atoms[count]),
                            float(element[0]),
                            float(element[1]),
                            float(element[2]),
                        )
                    )
                count += 1
            for element in linkatoms:
                # links are already in angstrom
                xyzfile.write(
                    "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                        str("H"), float(element[0]), float(element[1]), float(element[2])
                    )
                )
            #ofile.write("\n")
            #ofile.write("*" + "\n")


        with open("pointcharges.pc", 'w') as pfile:

            with open(pcffile) as ifile:
                c = 0
                for line in ifile:
                    match = re.search(
                        r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                    )
                    if match:
                        c += 1
            
            with open(pcffile) as ifile:
                pfile.write(str(int(c)) + "\n")
                for line in ifile:
                    match = re.search(
                        r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                    )
                    if match:       
                        pfile.write(
                            "{:>12.6f} {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                                float(match.group(4)),
                                float(match.group(1)),
                                float(match.group(2)),
                                float(match.group(3)),
                            )
                        )
        #ofile.write("\n")
        

    return serenityfile
'''

def make_orca_inp(qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    #qmmmtop = qmmmInputs.top                       SP activate next line
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra

    insert = ""
    oldinsert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step) ))
        if int(curr_step) > 1:
            oldinsert = str("." + str(int(curr_step) - 1))
    orcafile = str(jobname + insert + ".inp")
    gbwfile = str(jobname + insert + ".gbw")
    oldgbwfile = str(jobname + oldinsert + ".gbw")

    if nmaflag == 1:
        oldgbwfile = str(jobname + ".gbw")
    #chkfile = str(jobname + insert + ".chk")
    #oldchkfile = str(jobname + oldinsert + ".chk")

    #if nmaflag == 1:
    #    oldchkfile = str(jobname + ".chk")

    with open(orcafile, "w") as ofile:
        fullcoords = get_full_coords_angstrom(gro)
        atoms = get_atoms(qmmmtop, logfile)
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('top_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('QMMMTOPINFO')
            for i in atoms:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        with open('gro_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('GROINFO')
            for i in fullcoords:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        
        ofile.write("!" + str(method))
        if str(basis) != "NONE":
            ofile.write(" " + str(basis))
        ofile.write(" AutoAux ENGRAD PrintMOs Printbasis NoUseSym SCFCONV8 CHELPG" + "\n")
        if str(extra) != "NONE":
            ofile.write(str(extra) + "\n")
        if int(curr_step) == 0:
            ofile.write("!HUECKEL" + "\n")
        if int(curr_step) != 0 or nmaflag == 1:
            ofile.write("!MORead" + "\n")
            ofile.write("%moinp " + "\"" + oldgbwfile + "\"" +"\n")
        
        ofile.write("%maxcore" + " " + str(memory) + "\n")
        ofile.write("%pal" + "\n" + "nprocs" + " " + str(cores) + "\n" + "end" + "\n")
        

        #ofile.write("%moinp=" + gbwfile + "\n")
        #if int(curr_step) != 0 or nmaflag == 1:
        #    ofile.write("%OLDCHK=" + oldchkfile + "\n")
        #ofile.write("#P " + str(method))
        #if str(basis) != "NONE":
        #    ofile.write("/" + str(basis))
        #if str(extra) != "NONE":
        #    ofile.write(" " + str(extra))
        #if int(curr_step) != 0 or nmaflag == 1:
        #    ofile.write(" guess=read")
        #ofile.write(
        #    " nosymm gfinput gfprint force charge guess=huckel punch=derivatives iop(3/33=1,6/7=3) prop(field,read) pop=esp\n"
        #) #SP added 6/7=3 to print all MOs
        ofile.write("%pointcharges \"pointcharges.pc\" " + "\n\n")
        ofile.write(
            "\n#QM part of step "+str(curr_step)+"\n\n"
            + "*"
            + " "
            + "xyz"
            + " "
            + str(int(charge))
            + " "
            + str(int(multi))
            + "\n"
        )
        
        count = 0
        for element in fullcoords:
            if int(count + 1) in np.array(qmatomlist).astype(int):
                #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                ofile.write(
                    "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                        str(atoms[count]),
                        float(element[0]),
                        float(element[1]),
                        float(element[2]),
                    )
                )
            count += 1
        for element in linkatoms:
            # links are already in angstrom
            ofile.write(
                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                    str("H"), float(element[0]), float(element[1]), float(element[2])
                )
            )
        #ofile.write("\n")
        ofile.write("*" + "\n")


        with open("pointcharges.pc", 'w') as pfile:

            with open(pcffile) as ifile:
                c = 0
                for line in ifile:
                    match = re.search(
                        r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                    )
                    if match:
                        c += 1
            
            with open(pcffile) as ifile:
                pfile.write(str(int(c)) + "\n")
                for line in ifile:
                    match = re.search(
                        r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                    )
                    if match:       
                        pfile.write(
                            "{:>12.6f} {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                                float(match.group(4)),
                                float(match.group(1)),
                                float(match.group(2)),
                                float(match.group(3)),
                            )
                        )
        #ofile.write("\n")
        

    return orcafile



def make_g16_inp(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    #qmmmtop = qmmmInputs.top                       SP activate next line
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    pcffile = qmmmInputs.pcffile
    curr_step = qmmmInputs.qmmmparams.curr_step
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    nmaflag = qmmmInputs.nmaflag
    #inout = qmmmInputs.inout # Added by Nicola
    #if inout:# Added by Nicola
    #    inner = qmmmInputs.innerlist # Added by Nicola
    #    outer = qmmmInputs.outerlist # Added by Nicola
    #    ref = qmatomlist + inner + outer # Added by Nicola
    #    ref = list(set(ref)) # Added by Nicola
    #    ref.sort() # Added by Nicola
    #    memory_dict = reindexing.reindexing_memory(ref, outer)
    #    for i in range(len(qmatomlist)): # Added by Nicola
    #        qmatomlist[i] = memory_dict[qmatomlist[i]] # Added by Nicola
        #for i in range(len(linkatoms)):# Added by Nicola
        #    linkatoms[i] = memory_dict[linkatoms[i]]# Added by Nicola

    method = qmmmInputs.qmparams.method
    basis = qmmmInputs.qmparams.basis
    charge = qmmmInputs.qmparams.charge
    multi = qmmmInputs.qmparams.multi
    cores = qmmmInputs.qmparams.cores
    memory = qmmmInputs.qmparams.memory
    extra = qmmmInputs.qmparams.extra
    
    insert = ""
    oldinsert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step) ))
        if int(curr_step) > 1:
            oldinsert = str("." + str(int(curr_step) - 1))
    gaufile = str(jobname + insert + ".gjf")
    chkfile = str(jobname + insert + ".chk")
    oldchkfile = str(jobname + oldinsert + ".chk")

    if nmaflag == 1:
        oldchkfile = str(jobname + ".chk")

    with open(gaufile, "w") as ofile:
        fullcoords = get_full_coords_angstrom(gro)
        atoms = get_atoms(qmmmtop, logfile)
        original_stdout = sys.stdout # Save a reference to the original standard output
        with open('top_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('QMMMTOPINFO')
            for i in atoms:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        with open('gro_info.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('GROINFO')
            for i in fullcoords:
                print(str(i))
            sys.stdout = original_stdout # Reset the standard output to its original value
        ofile.write("%NPROCSHARED=" + str(cores) + "\n")
        ofile.write("%MEM=" + str(memory) + "MB\n")
        ofile.write("%CHK=" + chkfile + "\n")
        if int(curr_step) != 0 or nmaflag == 1:
            ofile.write("%OLDCHK=" + oldchkfile + "\n")
        ofile.write("#P " + str(method))
        if str(basis) != "NONE":
            ofile.write("/" + str(basis))
        if str(extra) != "NONE":
            ofile.write(" " + str(extra))
        if int(curr_step) != 0 or nmaflag == 1:
            ofile.write(" guess=read")
        ofile.write(
            " nosymm gfinput gfprint force charge guess=huckel punch=derivatives iop(3/33=1,6/7=3) prop(field,read) pop=esp\n"
        ) #SP added 6/7=3 to print all MOs
        ofile.write(
            "\nQM part of step "+str(curr_step)+"\n\n"
            + str(int(charge))
            + " "
            + str(int(multi))
            + "\n"
        )
        count = 0
        for element in fullcoords:
            if int(count + 1) in np.array(qmatomlist).astype(int):
                #print(str(count)+" "+str(element)+" "+str(len(atoms)))
                ofile.write(
                    "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                        str(atoms[count]),
                        float(element[0]),
                        float(element[1]),
                        float(element[2]),
                    )
                )
            count += 1
        for element in linkatoms:
            # links are already in angstrom
            ofile.write(
                "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                    str("H"), float(element[0]), float(element[1]), float(element[2])
                )
            )
        ofile.write("\n")

        with open(pcffile) as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                )
                if match:
                    ofile.write(
                        "{:>12.6f} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                            float(match.group(4)),
                        )
                    )
        ofile.write("\n")

        with open(pcffile) as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                )
                if match:
                    ofile.write(
                        "{:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            float(match.group(1)),
                            float(match.group(2)),
                            float(match.group(3)),
                        )
                    )
        ofile.write("\n\n")

    return gaufile

def make_gmx_inp(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    gro = qmmmInputs.gro
    qmmmtop = qmmmInputs.qmmmtop
    qmatomlist = qmmmInputs.qmatomlist
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    if qmmmInputs.inout:
        inout = True
    else:
        inout= False
    
    if qmmmInputs.mmparams.gmxplus == 0:
        gmxplus = False
    elif qmmmInputs.mmparams.gmxplus == 1:
        gmxplus = True

    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))
    mdpname = str(jobname + ".mdp")
    groname = str(jobname + ".boxlarge.g96")
    ndxname = str(qmmmtop + ".ndx")
    tprname = str(jobname + insert + ".tpr")
    rcoulomb = qmmmInputs.mmparams.rcoulomb
    rvdw = qmmmInputs.mmparams.rvdw
    nbradius = get_nbradius(gro)
    write_mdp(mdpname, rcoulomb, rvdw, nbradius, inout, gmxplus)
    update_gro_box(gro, groname, nbradius, logfile)
    subprocess.call(
        [
            prefix,
            "grompp",
            "-p",
            str(qmmmtop),
            "-c",
            str(groname),
            "-n",
            str(ndxname),
            "-f",
            str(mdpname),
            "-o",
            str(tprname),
            "-backup",
            "no",
        ]
    )
    subprocess.call(["rm", "mdout.mdp"])
    return tprname

def write_mdp(mdpname, rcoulomb, rvdw, nbradius, inout, gmxplus):#gmxplus = True added by Nicola
    if rcoulomb == 0:
        rcoulomb = nbradius
    if rvdw == 0:
        rvdw = rcoulomb
    with open(mdpname, "w") as ofile:
        ofile.write(
            "title               =  Yo\ncpp                 =  /usr/bin/cpp\nconstraints         =  none\nintegrator          =  md\ndt                  =  0.001 ; ps !\nnsteps              =  1\nnstcomm             =  0\nnstxout             =  1\nnstvout             =  1\nnstfout             =   1\nnstlog              =  1\nnstenergy           =  1\nnstlist             =  1\nns_type             =  grid\nrlist               =  "
        )
        ofile.write(str(float(rcoulomb)))
        if not gmxplus:#Added by Nicola
            ofile.write(
                "\ncutoff-scheme = group\ncoulombtype    =  cut-off\nrcoulomb            =  "
            )
        else:#Added by Nicola
            ofile.write(
                "\ncutoff-scheme = verlet\ncoulombtype    =  cut-off\nrcoulomb            =  "
            )
        ofile.write(str(float(rcoulomb)))
        ofile.write("\nrvdw                =  ")
        ofile.write(str(float(rvdw)))
        '''if not gmxplus: #Added by Nicola
            if inout:
                ofile.write(
                "\nTcoupl              =  no\nfreezegrps          =  OUTER\nfreezedim           =  Y Y Y\nenergygrps          =  QM INNER OUTER\nenergygrp-excl = QM QM INNER OUTER OUTER OUTER\nPcoupl              =  no\ngen_vel             =  no\n"
                )
            else:
                ofile.write(
                "\nTcoupl              =  no\nenergygrps          =  QM\nenergygrp-excl = QM QM\nPcoupl              =  no\ngen_vel             =  no\n"
                )'''
        #else:#Added by Nicola
        #    if inout:
        ofile.write(
        "\nTcoupl              =  no\nPcoupl              =  no\ngen_vel             =  no\n"
        )
        #    else:
        #        ofile.write(
        #        "\nTcoupl              =  no\nPcoupl              =  no\ngen_vel             =  no\n"
        #        )
            #if inout:
            #    ofile.write(
            #    "\nTcoupl              =  no\nfreezegrps          =  OUTER\nfreezedim           =  Y Y Y\nPcoupl              =  no\ngen_vel             =  no\n"
            #    )
            #else:
            #    ofile.write(
            #    "\nTcoupl              =  no\nPcoupl              =  no\ngen_vel             =  no\n"
            #    )

def update_gro_box(gro, groname, nbradius, logfile):
    with open(groname, "w") as ofile:
        with open(gro) as ifile:
            logger(
                logfile,
                str(
                    "Finding a larger .gro box size to avoid problems with .mdp input...\n"
                ),
            )
            for line in ifile:
                ofile.write(line)
                match = re.search(r"^BOX\s*\n", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\d*\.\d+)\s+(\d*\.\d+)\s+(\d*\.\d+)",
                            line,
                            flags=re.MULTILINE,
                        )
                        if not match:
                            logger(
                                logfile,
                                "\n\nError: In "
                                + str(gro)
                                + " box vectors were expected but not found. Exiting. Line was:\n",
                            )
                            logger(logfile, line)
                            exit(1)
                        else:
                            bv = [
                                float(match.group(1)) + 10.0 * nbradius,
                                float(match.group(2)) + 10.0 * nbradius,
                                float(match.group(3)) + 10.0 * nbradius,
                            ]
                            ofile.write(
                                " {:>15.9f} {:>15.9f} {:>15.9f}\nEND\n".format(
                                    float(bv[0]), float(bv[1]), float(bv[2])
                                )
                            )
                        break
                    break
    logger(logfile, str("Done.\n"))

#Propagated
#initstep = stepsize
def propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile, scan_flag=False, scan_atoms=[]):
    dispvec = []
    maxforce = 0.0
    clean_force = make_clean_force(total_force)
    old_clean_force = make_clean_force(last_forces)

    if scan_flag :
        scan_type = len(scan_atoms)
        if scan_type == 2 :
            clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'R')
            if (propagator == "BFGS") and (len(old_clean_force) > 1):
                old_clean_force = scan_adj_force(old_clean_force, xyzq, scan_atoms, 'R')

        elif scan_type == 3 :
            clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'A')

        elif scan_type == 4 :
            clean_force = scan_adj_force(clean_force, xyzq, scan_atoms, 'B')

    maxatom = -1
    maxcoord = -1
    check_force = list(_flatten(clean_force))
    #search index of max info
    for i in range(0, int(len(check_force) / 3)):
        for j in range(0, 3):
                if abs(float(check_force[i * 3 + j])) > abs(maxforce):
                    maxforce = float(check_force[i * 3 + j])
                    maxatom = i
                    maxcoord = j
    logger(
        logfile,
        str(
            "Maximum force is "
            + str(float(maxforce))
            + " a.u. at coord "
            + str(int(maxatom) + 1)
            + "/"
            + str(int(maxcoord) + 1)
            + ".\n"
        ),
    )
    #All use steep for first propagation
    if curr_step <= 1:
        logger(logfile, "Propagate with steepest descent for first step...\n")
        for element in clean_force:
            dispvec.append(
                [
                    float(element[0]) * float(stepsize) / abs(float(maxforce)),
                    float(element[1]) * float(stepsize) / abs(float(maxforce)),
                    float(element[2]) * float(stepsize) / abs(float(maxforce)),
                ]
            )
    else:
        if propagator == "STEEP":
            logger(logfile, "Propagate with steepest descent...\n")
            for element in clean_force:
                dispvec.append(
                    [
                        float(element[0]) * float(stepsize) / abs(float(maxforce)),
                        float(element[1]) * float(stepsize) / abs(float(maxforce)),
                        float(element[2]) * float(stepsize) / abs(float(maxforce)),
                    ]
                )
            corr_length = np.array(total_force)

        elif propagator == "CONJGRAD" :
            logger(logfile, "Propagate with conjugate gradient...\n")
            # Fletcher-Reeves
            _flattened = list(_flatten(clean_force))
            corr_fac = np.array(_flattened).dot(np.array(_flattened))
            _flattened = list(_flatten(old_clean_force))
            corr_fac /= np.array(_flattened).dot(np.array(_flattened))

            #2222
            corr_length = np.array(total_force)
            if curr_step > 1:   #0 step for the inital SP 
                corr_length = np.array(total_force) + corr_fac * corr_length

            for i in range(len(total_force)):
                dispvec.append(
                    [
                        float(stepsize) * corr_length[i][0],
                        float(stepsize) * corr_length[i][1],
                        float(stepsize) * corr_length[i][2],
                    ]
                )

            logger(
                logfile,
                str(
                    "Effective step at maximum force coord is "
                    + str(float(dispvec[maxatom][maxcoord]))
                    + " a.u.\n"
                ),
            )
        elif propagator == "BFGS":
            logger(logfile, "Propagate with BFGS method...\n")
            coords = np.array(new_xyzq)[:, 0:3]
            old_coords = np.array(xyzq)[:, 0:3]
            old_hessian = np.loadtxt("bfgs_hessian.txt")
            gradient = np.array(total_force)
            old_gradient = np.array(last_forces)
            hessian, hesseig, warning_flag = get_approx_hessian(
                coords, old_coords, gradient, old_gradient, old_hessian, logfile
            )
            np.savetxt("bfgs_hessian.txt", hessian)

            coords = coords.reshape(
                3 * len(coords), 1
            )  # reshape coords to use in dot products
            gradient = gradient.reshape(
                3 * len(gradient), 1
            )  # reshape grad to use in dot products

            direc = -np.linalg.inv(hessian).dot(gradient)
            direc = direc.reshape(int(len(coords) / 3), 3)

            for i in range(len(total_force)):
                dispvec.append(
                    [
                        float(stepsize) * direc[i][0],
                        float(stepsize) * direc[i][1],
                        float(stepsize) * direc[i][2],
                    ]
                )
    

    return dispvec

def g96_to_gro(inp,out,logfile):
    import re
    with open(out,"w") as ofile:
        n_a=0
        with open(inp) as ifile:
            for line in ifile:
                match=re.search(r'^POSITION', line,flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match=re.search(r'^END', line,flags=re.MULTILINE)
                        if match:
                           break
                        n_a+=1
                    break
        ofile.write("GRoups of Organic Molecules in ACtion for Science\n")
        if n_a>99999:
            ofile.write(str(n_a)+"\n")
        else:
            ofile.write("{:>5d}\n".format(n_a))
        with open(inp) as ifile:
            for line in ifile:
                match=re.search(r'^POSITION', line,flags=re.MULTILINE)
                if match:
                   for line in ifile:
                       match=re.search(r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)", line,flags=re.MULTILINE)
                       if match:
                           resid=int(match.group(1))
                           atomid=int(match.group(4))
                           while resid>99999:
                                   resid-=100000
                           while atomid>99999:
                                   atomid-=100000
                           ofile.write("{:>5d}{:<5}{:>5}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(resid,match.group(2),match.group(3),atomid,float(match.group(5)),float(match.group(6)),float(match.group(7))))
                           continue
                       match=re.search(r"^END", line,flags=re.MULTILINE)
                       if not match:
                           logger(logfile,str("Unexpected entry in .g96 file. Exiting. Last line:\n"))
                           logger(logfile,line)
                           exit(1)
                       else:
                           break
                match=re.search(r"^BOX", line,flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^\s*(\d*\.\d+)\s+(\d*\.\d+)\s+(\d*\.\d+)",line,flags=re.MULTILINE)
                        if not match:
                            logger(logfile,str("Unexpected entry in .g96 file. Exiting. Last line:\n"))
                            logger(logfile,line)
                            exit(1)
                        else:
                            ofile.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(float(match.group(1)),float(match.group(2)),float(match.group(3))))
                            break
                    break

def make_g96_inp(dispvec, gro, new_gro, logfile):

    with open(new_gro, "w") as ofile:
        with open(gro) as ifile:
            counter = 0
            for line in ifile:
                ofile.write(line)
                counter += 1
                if counter == 4:
                    break
            counter = 0
            for line in ifile:
                match = re.search(
                    r"^(.{5})\s(.{5})\s(.{5})\s(.{6})\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)\s*([-]*\d+\.*\d*)",
                    line,
                    flags=re.MULTILINE,
                )
                if not match:
                    ofile.write(line)
                    logger(
                        logfile,
                        str(
                            "Successfully wrote "
                            + str(int(counter))
                            + " atoms to new g96 file.\n"
                        ),
                    )
                    break
                else:
                    dispx = dispvec[counter][0] * 0.052917721
                    dispy = dispvec[counter][1] * 0.052917721
                    dispz = dispvec[counter][2] * 0.052917721
                    ofile.write(
                        str(match.group(1))
                        + " "
                        + str(match.group(2))
                        + " "
                        + str(match.group(3))
                        + " "
                        + str(match.group(4))
                        + " {:>15.9f} {:>15.9f} {:>15.9f}\n".format(
                            float(match.group(5)) + float(dispx),
                            float(match.group(6)) + float(dispy),
                            float(match.group(7)) + float(dispz),
                        )
                    )
                    counter += 1
            for line in ifile:
                ofile.write(line) 
    
    ifile.close()
    ofile.close()
    logger(logfile, str("Made new coordinates in file " + str(new_gro) + ".\n"))   

def qmmm_prep(qmmmInputs):
    gro = qmmmInputs.gro
    top = qmmmInputs.top
    jobname = qmmmInputs.qmmmparams.jobname
    qmatomlist = qmmmInputs.qmatomlist
    connlist = qmmmInputs.connlist
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    gmxtop = qmmmInputs.pathparams.gmxtop 
    charge = qmmmInputs.qmparams.charge
    geo = make_pcf.readg96(gro)
    curr_step = qmmmInputs.qmmmparams.curr_step

    logger(logfile, "List of molecules...\n")
    # 6277
    mollist = make_pcf.readmols(top)
    logger(logfile, "Done.\n")
    logger(logfile, "Reading charges...\n")
    chargevec = []
    for element in mollist:
        chargevec.extend(make_pcf.readcharges(element, top, gmxtop))
    logger(logfile, "Done.\n")
    #if qmmmInputs.inout and qmmmInputs.mmparams.gmxplus == 0:
    #    new_xyzq = make_xyzq_io(geo, chargevec, qmmmInputs.outerlist)
    #else:
    new_xyzq = make_xyzq(geo, chargevec)
    logger(logfile, str("Made new xyzq matrix.\n"))
    logger(
        logfile,
        "Preparing the point charge field for a numerically optimized charge shift...\n",
    )
    (
        qmcoordlist,
        m1list,
        m2list,
        updated_chargelist,
    ) = prep_pcf.prepare_pcf_for_shift_fieldsonly(
        new_xyzq, qmatomlist, charge, connlist
    )
    logger(logfile, "Done.\n")
    new_links = get_linkatoms_ang(
        new_xyzq, qmatomlist, m1list, connlist, linkatoms
    )
    logger(logfile, str("Updated positions of link atoms.\n"))
    filename = jobname
    if curr_step > 0:
        filename += "." + str(int(curr_step))
    if not os.path.isfile(str(filename + ".pointcharges")):
        logger(logfile, "Shifting...\n")
        final_pcf.generate_charge_shift_fieldsonly(
            updated_chargelist, m1list, qmcoordlist, m2list, filename, basedir
        )
        logger(logfile, str("Made new PCF file.\n"))
    else:
        logger(
            logfile,
            "NOTE: Shifting omitted due to "
            + str(filename + ".pointcharges")
            + " being an existing file!\n",
        )
    logger(logfile, "Done.\n")

    #Update new output
    qmmmInputs.xyzq = new_xyzq
    qmmmInputs.m1list = m1list
    qmmmInputs.m2list = m2list
    qmmmInputs.linkatoms = new_links


    return qmmmInputs

#Database
def databasecorrection(energy_or_force, cut, dist, qmmmInputs):
    
    forcefield = qmmmInputs.mmparams.fField
    method = qmmmInputs.qmparams.method
    basisset = qmmmInputs.qmparams.basis
    fit = qmmmInputs.qmmmparams.databasefit
    basedir = qmmmInputs.basedir
    logfile = qmmmInputs.logfile

    conn = sqlite3.connect(basedir + "/correction_database/database.sqlite")

    # check if method exist in database
    method_set = conn.cursor()
    method_set.execute(
        "SELECT * FROM "
        + cut
        + ' WHERE forcefield="'
        + forcefield
        + '" AND method="'
        + method
        + '" AND basisset ="'
        + basisset
        + '"'
    )
    method_set_value = method_set.fetchall()
    if len(method_set_value) == 0:
        cut = "aminoacid_CACB"
        forcefield = "amberGS"
        method = "CAM-B3LYP"
        basisset = "6-31++G**"
        logger(
            logfile,
            "Unexisted method in correction database, changing to default correction method...\n",
        )

    c = conn.cursor()
    c.execute(
        "SELECT * FROM "
        + cut
        + ' WHERE forcefield="'
        + forcefield
        + '" AND method="'
        + method
        + '" AND basisset ="'
        + basisset
        + '"'
    )
    db_values = c.fetchall()[0]

    conn.close()
    returnvalue = 0
    if len(db_values) > 0:
        if energy_or_force == "ENERGY":
            if fit == "POLY":
                returnvalue = (
                    db_values[5] * dist * dist * dist
                    + db_values[6] * dist * dist
                    + db_values[7] * dist
                    + db_values[8]
                )
            elif fit == "MORSE":
                returnvalue = (
                    db_values[9]
                    * (
                        np.exp(-2 * db_values[10] * (dist - db_values[11]))
                        - 2 * np.exp(-db_values[10] * (dist - db_values[11]))
                    )
                    + db_values[12]
                )
            elif fit == "NO":
                returnvalue = 0
                logger(logfile, "No energy correction.\n")
        elif energy_or_force == "FORCES":
            if fit == "POLY":
                returnvalue = (
                    db_values[13] * dist * dist * dist
                    + db_values[14] * dist * dist
                    + db_values[15] * dist
                    + db_values[16]
                )
            elif fit == "MORSE":
                returnvalue = (
                    db_values[17]
                    * (
                        np.exp(-2 * db_values[18] * (dist - db_values[19]))
                        - 2 * np.exp(-db_values[18] * (dist - db_values[19]))
                    )
                    + db_values[20]
                )
            elif fit == "NO":
                returnvalue = 0
                logger(logfile, "No force correction.\n")
    return returnvalue
'''
def check_error_in_orcalog(log_filename):
    error_string = "ORCA finished by error termination in"
    mem_string1 = "Not enough memory"
    mem_string2 = "not enough memory"
    mem_string3 = "Out of memory"
    mem_string4 = "out of memory"
    mem_found = False
    error_found = False
    with open(log_filename, "r") as log_file:
        lines = log_file.readlines()
        lines.reverse()

        for line in lines:
            if mem_string1 in line or mem_string2 in line or mem_string3 in line or mem_string4 in line:
                mem_found = True
                break
        
        for line in lines:
            if error_string in line:
                error_found = True
                if mem_found:
                    return False, error_found, mem_found 
                else:
                    return True, error_found, mem_found
    return False, error_found, mem_found
'''
def check_error_in_orcalog(log_filename):
    error_string = "ORCA finished by error termination in"
    mem_strings = [
        "Not enough memory", 
        "not enough memory", 
        "Out of memory", 
        "out of memory"
    ]
    mem_found = False
    error_found = False

    with open(log_filename, "r") as log_file:
        lines = log_file.readlines()
        lines.reverse()

        for line in lines:
            if any(mem_string in line for mem_string in mem_strings):
                mem_found = True
            if error_string in line:
                error_found = True
                break

    return error_found, mem_found

def check_error_in_serenitylog(log_filename):
    with open(log_filename, "r") as log_file:
        lines = log_file.readlines()
        lines.reverse()
        for line in lines:
            if 'Serenity Crashed' in line:
                return True
            else:
                return False

#Run program command
def run_TM(qmfile, qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    #2222
    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))
    
    log_name = str(jobname + insert + ".inp.log")
    #engrad_name = str(jobname + insert + ".engrad")
    if not os.path.isfile(str(qmfile) + ".log"):
        logger(logfile, "Running turbomole file.\n")
        command = f" ridft > {qmfile[:-3]}log && rdgrad >> {qmfile[:-3]}log"
        subprocess.call(command, shell=True)
        #subprocess.call([orcacmd, str(qmfile)])
        logname = qmfile[:-3]
        logname += "log"
        subprocess.call(["mv", logname, log_name])
        subprocess.call(["mv", "control", qmfile])
        subprocess.call(f"cat energy gradient {qmfile} >> {log_name}",shell=True)
        subprocess.call(["mv", "basis", str(qmfile[:-3] + "basis")])
        subprocess.call(["mv", "auxbasis", str(qmfile[:-3] + "auxbasis")])
        subprocess.call(["mv", "coord", str(qmfile[:-3] + "coord")])
        subprocess.call(["mv", "*mos", str(qmfile[:-3] + "mos")])
        subprocess.call(["mv", "str", str(qmfile[:-3] + "str")])
        subprocess.call("rm energy gradient", shell=True)
        if os.path.isfile(log_name):
            logger(logfile, "Turbomole done.\n")
        else:
            logger(logfile, "Turbomole job was not done. Exiting.\n")
            exit(1)
    else:
        logger(logfile, "NOTE: Using existing serenity files, skipping calculation for this step.\n")


def run_serenity(qmfile, qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    serenitycmd = qmmmInputs.pathparams.serenitycmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    #2222
    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))
    
    log_name = str(jobname + insert + ".inp.log")
    #engrad_name = str(jobname + insert + ".engrad")
    if not os.path.isfile(str(qmfile) + ".log"):
        logger(logfile, "Running serenity file.\n")
        command = f"{serenitycmd} {qmfile} > {qmfile[:-3]}log"
        subprocess.call(command, shell=True)
        #subprocess.call([orcacmd, str(qmfile)])
        logname = qmfile[:-3]
        logname += "log"
        subprocess.call(["mv", logname, log_name])
        if os.path.isfile(log_name):
            error_found = check_error_in_serenitylog(log_name)
            if error_found:
                logger(logfile, "ERROR was detected, serenity run was aborted, Exiting serenity.\n")
                exit(1)
            else:
                logger(logfile, "Serenity done.\n")
    else:
        logger(logfile, "NOTE: Using existing serenity files, skipping calculation for this step.\n")

        

def run_orca(qmfile, qmmmInputs): # Added by Nicola
    jobname = qmmmInputs.qmmmparams.jobname
    orcacmd = qmmmInputs.pathparams.orcacmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    #2222
    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))
    
    log_name = str(jobname + insert + ".inp.log")
    engrad_name = str(jobname + insert + ".engrad")
    while True:
        if not os.path.isfile(str(qmfile) + ".log"):
            logger(logfile, "Running orca file.\n")
            command = f"{orcacmd} {qmfile} > {qmfile[:-3]}log"
            subprocess.call(command, shell=True)
            #subprocess.call([orcacmd, str(qmfile)])
            logname = qmfile[:-3]
            logname += "log"
            subprocess.call(["mv", logname, log_name])
            if os.path.isfile(log_name):
                error_found, mem_found = check_error_in_orcalog(log_name)
                if error_found:
                    if mem_found:
                        logger(logfile, "ERROR: Detected error termination in ORCA log due to insufficient memory. Exiting ORCA.\n")
                        exit(1)
                    else:
                        logger(logfile, "ERROR: Detected error termination in ORCA log. Rerunning ORCA.\n")
                        subprocess.call(f"rm *.tmp {log_name}", shell=True)
                        continue
                else:
                    subprocess.call(["mv", f"{qmfile[:-4]}.engrad", engrad_name])
                    logger(logfile, "ORCA Done.\n")
                    break
        else:
            logger(logfile, "NOTE: Using existing ORCA files, skipping calculation for this step.\n")
            break

    if not os.path.isfile(engrad_name):
        if not os.path.isfile(qmfile[:-4] + ".engrad"):
            logger(
                logfile,
                "No engrad file was created by the last Orca run! Exiting.\n",
            )
            exit(1)
        subprocess.call(["mv", qmfile[:-4] + ".engrad", engrad_name])
        logger(
            logfile,
            "WARNING: Had to rename engrad file but not the log file. MAKE SURE THAT THE ENGRAD FILE FITS TO THE LOG FILE!\n",
        )



def run_g16(qmfile, qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    g16cmd = qmmmInputs.pathparams.g16cmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    #2222
    insert = ""
    if int(curr_step) > 0:
        insert = str("." + str(int(curr_step)))


    if not os.path.isfile(str(qmfile) + ".log"):
        logger(logfile, "Running G16 file.\n")
        subprocess.call([g16cmd, str(qmfile)])
        logname = qmfile[:-3]
        logname += "log"
        subprocess.call(["mv", logname, str(jobname + insert + ".gjf.log")])
        subprocess.call(["mv", "fort.7", str(jobname + insert + ".fort.7")])
        logger(logfile, "G16 Done.\n")
    else:
        logger(
            logfile,
            "NOTE: Using existing G16 files, skipping calculation for this step.\n",
        )
    if not os.path.isfile(jobname + insert + ".fort.7"):
        if not os.path.isfile("fort.7"):
            logger(
                logfile,
                "No fort.7 file was created by the last Gaussian run! Exiting.\n",
            )
            exit(1)
        subprocess.call(["mv", "fort.7", str(jobname + insert + ".fort.7")])
        logger(
            logfile,
            "WARNING: Had to rename fort.7 file but not the log file. MAKE SURE THAT THE FORT.7 FILE FITS TO THE LOG FILE!\n",
        )

def run_gmx(mmfile, qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile

    insert = ""
    if int(curr_step) != 0:
        insert = str("." + str(int(curr_step)))

    logger(logfile, "Running Gromacs file.\n")
    trrname = str(jobname + insert + ".trr")
    xtcname = str(jobname + insert + ".xtc")
    outname = str(jobname + insert + ".out.gro")
    gmxlogname = str(jobname + insert + ".gmx.log")
    edrname = str(jobname + insert + ".edr")
    subprocess.call(
        [
            prefix,
            "mdrun",
            "-s",
            mmfile,
            "-o",
            trrname,
            "-c",
            outname,
            "-x",
            xtcname,
            "-g",
            gmxlogname,
            "-e",
            edrname,
            "-backup",
            "no",
        ]
    )
    subprocess.call(["rm", outname])

    return edrname

# Write output file
    #old output
def old_write_output(energies, total_force, curr_step ):
    qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    energy_file = "oenergy.txt"
    
    if curr_step == 0:
        file_flag = "w"
    else:
        
        file_flag = "a+"
    
    oenergy = open(energy_file, file_flag)
    oenergy.write("Step: %d\n" % curr_step)
    oenergy.write("QM energy: %.8f\n" % qmenergy)
    oenergy.write("MM energy: %.8f\n" % mmenergy)
    oenergy.write("Link energy: %.8f\n" % linkcorrenergy)
    oenergy.write("Total energy: %.8f\n"% total_energy)
    oenergy.write("--------------------------------------------\n")
    oenergy.close()

    oforce = open("oforce.txt", file_flag)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write('\n')
    oforce.close()
    #test output

def write_test(qmmmInputs, curr_step):
    qmenergy, mmenergy, linkcorrenergy, total_energy = qmmmInputs.energies
    total_force = qmmmInputs.forces
    #qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    energy_file = "testE.txt"
    
    if curr_step == 0:
        file_flag = "w"
    else:
        
        file_flag = "a+"
    
    oenergy = open(energy_file, file_flag)
    oenergy.write("Step: %d\n" % curr_step)
    oenergy.write("QM energy: %.8f\n" % qmenergy)
    oenergy.write("MM energy: %.8f\n" % mmenergy)
    oenergy.write("Link energy: %.8f\n" % linkcorrenergy)
    oenergy.write("Total energy: %.8f\n"% total_energy)
    oenergy.write("--------------------------------------------\n")
    oenergy.close()

    oforce = open("testF.txt", file_flag)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write('\n')
    oforce.close()

def write_output(energies, total_force, curr_step, energy_file="oenergy.txt", forces_file="oforces.txt"):
    qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    #energy_file = "oenergy.txt"
    
    if curr_step == 0:
        file_flag = "w"
    else:
        file_flag = "a+"

    oenergy = open(energy_file, file_flag)
    
    if curr_step == 0:
        oenergy.write("Step\tQM\t\tMM\t\tLink\t\tTotal\n")
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (curr_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    else:
        
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (curr_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    oenergy.close()

    oforce = open(forces_file, file_flag)
    oforce.write("Step%d\n"%curr_step)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write("\n")
    oforce.close()

def write_statechar(char, curr_step):
    state_file = "ostatechar.txt"

    if curr_step == 0:
        file_flag = "w"
    else:
        file_flag = "a+"

    out = open(state_file, file_flag)

    if curr_step == 0:
        out.write("Statecharacter\n")
        for item in char:
            out.write("%d\t%d\t%f\n" % (item[0],item[1],item[2]))
    else:
        for item in char:
            out.write("%d\t%d\t%f\n" % (item[0],item[1],item[2]))
    out.close()

def write_dispvec(dispvec, curr_step, count_trash):
    dispvec_file = "dispvec_"+str(curr_step)+"_"+str(count_trash)+".txt"
    file_flag = "w"

    out = open(dispvec_file, file_flag)
    out.write("current step"+str(curr_step)+"rejection time"+str(count_trash)+"\n")
    for item in dispvec:
        out.write("%f\t%f\t%f\n" % (item[0],item[1],item[2]))
    out.close()

def write_scan(energies, total_force, scan_step, energy_file='oenergy.txt', forces_file='oforces.txt'):#SP
    qmenergy, mmenergy, linkcorrenergy, total_energy = energies
    if scan_step == 1:
        file_flag = "w"
    else:
        file_flag = "a+"

    oenergy = open(energy_file, file_flag)

    if scan_step == 1:
        oenergy.write("Scan_step\tQM\t\tMM\t\tLink\t\tTotal\n")
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (scan_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    else:
        oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (scan_step, qmenergy, mmenergy, linkcorrenergy, total_energy))
    oenergy.close()

    #forces_file = "oforces.txt"
    oforce = open(forces_file, file_flag)
    oforce.write("Scan_step%d\n"%scan_step)
    for i in range(len(total_force)):
        oforce.write(
            "%d %.8f %.8f %.8f\n"
            % ((i + 1), total_force[i][0], total_force[i][1], total_force[i][2])
        )
    oforce.write("\n")
    oforce.close()

# Read output file
def read_oenergy(step, filename = 'oenergy.txt'):
    oenergy = np.loadtxt(filename,dtype=float,skiprows=1,usecols=(1,2,3,4))
    if  isinstance(oenergy[step], np.ndarray):        # if one opt step is enough, oenergy[step] is just the last value AJ
        return oenergy[step]
    else:
        return oenergy

def read_oforces(step, filename = 'oforces.txt'):
    file_info = open(filename, 'r')
    step_index = []
    for i, line in enumerate(file_info.readlines()): 
        if 'Step' in line:
            step_index.append(i)
    steps = len(step_index)
    with open(filename,'r') as file_object:
        line = file_object.readlines()
    for i in range(len(step_index)):
        if i < len(step_index)-1:
            outlines = step_index[i+1] - step_index[i]
        else:
            outlines = len(line) - step_index[i] - 2                #why -2 in old version ??
    forces = []
    for i in range(outlines):
        force = line[step_index[step]+i+1].split()[1:]
        for j in range(3):
            force[j] = float(force[j])
        forces.append(force)
    return np.array(forces)

# Archive & remove
def archive(jobname, curr_step):
    if curr_step == 0 :
        insert = '' 
    else:
        insert = str("." + str(curr_step))
        output_filename = str(jobname + insert + ".tar")
        archive_files = str(jobname + insert + "*")
        exclude_files = "*.tar"
        subprocess.call("tar --exclude=%s -czf %s %s --remove-files "%(exclude_files,output_filename, archive_files), shell=True)

def archive_trash(jobname, curr_step, cycle_step):                #Simon, if you want to check the counterproductive steps
    if curr_step == 0 :
        insert = ''
    else:
        insert = str("." + str(curr_step))
    output_filename = str(jobname + insert +"."+ str(cycle_step) + ".tar")
    archive_files = str(jobname + insert + "*")
    exclude_files = "*.tar"
    subprocess.call("tar --exclude=%s -czf %s %s --remove-files "%(exclude_files,output_filename, archive_files), shell=True)

def archive_first_step(jobname):
    output_filename = str("step_zero.tar")
    archive_files = str(jobname+"*")
    subprocess.call("tar -czf %s %s"%(output_filename, archive_files), shell=True)

def remove_files(jobname, curr_step, tar=False):
    if curr_step == 0 :
        insert = '' 
    else:
        insert = str("." + str(curr_step))

    extend = np.array([".g96",".boxlarge.g96",".chk",".edr",".edr.xvg",".fort.7",".gjf",".gjf.log",
        ".gmx.log", ".mdp", ".pointcharges", ".qmmm.top",".qmmm.top.ndx", ".tpr", ".trr", ".xvg"])

    for i in range(len(extend)):
        subprocess.call("rm -f %s"%str(jobname + insert + extend[i]), shell=True)

    if tar:
        subprocess.call("rm -f %s"%str(jobname + insert + '.tar'), shell=True)

#Energy
def get_qmenergy(qmfile, qmmmInputs):
    qmprog = qmmmInputs.qmparams.program
    extra_string = qmmmInputs.qmparams.extra
    pcffile = qmmmInputs.pcffile
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    if qmprog == "SERENITY":
        fde = qmmmInputs.qmparams.fde
        tddft = qmmmInputs.qmparams.tddft

    logger(logfile, "Extracting QM energy.\n")
    qmenergy = 0.0
    qm_corrdata = []
    if str(qmprog) == "G16":
        with open(str(qmfile + ".log")) as ifile:
            for line in ifile:
                match = []
                match2 = []
                match2 = re.search(
                    r"\sTD[=(\s]", extra_string.upper(), flags=re.MULTILINE
                )
                if not match2:
                    match2 = re.search(
                        r"^TD[=(\s]", extra_string.upper(), flags=re.MULTILINE
                    )
                if not match2:
                    match2 = re.search(
                        r"\sTD$", extra_string.upper(), flags=re.MULTILINE
                    )
                if not match2:
                    match2 = re.search(
                        r"^TD$", extra_string.upper(), flags=re.MULTILINE
                    )
                if not match2:
                    match = re.search(
                        r"^\s*SCF\s*Done:\s*E\(\S+\)\s*\=\s*([-]*\d+\.\d+)",
                        line,
                        flags=re.MULTILINE,
                    )
                else:
                    match = re.search(
                        r"^\s*Total\s*Energy,\s*E\(\S+\)\s*\=\s*([-]*\d+\.\d+)",
                        line,
                        flags=re.MULTILINE,
                    )
                if match:
                    logger(logfile, "Obtaining charge self-interaction...\n")
                    pcf_self_pot = read_pcf_self(qmfile)
                    logger(
                        logfile, "Done: {:>20.10f} a.u.\n".format(float(pcf_self_pot))
                    )
                    # G16 energy needs to be corrected for self potential of PCF
                    qmenergy = float(match.group(1)) - float(pcf_self_pot)
                match = re.search(r"^\s*ESP\s*charges:", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        break
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\d+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            qm_corrdata.append(
                                [
                                    int(match.group(1)),
                                    match.group(2),
                                    float(match.group(3)),
                                ]
                            )
                        else:
                            break
                    break
    if str(qmprog) == "ORCA": # Added by Nicola
        with open(str(qmfile[:-4] + ".engrad")) as ifile:
            lines = ifile.readlines()
            for line in range(len(lines)):
                if '# The current total energy in Eh' in lines[line]:
                    qmenergy = float(lines[line + 2].strip())# - float(calculate_nuclear_self_interaction(qmfile))
                    break
                else:
                    continue
        '''
        #Alternative
        with open(str(qmfile + ".log")) as ifile:
            lines = ifile.readlines()
            for line in lines:
                if 'FINAL SINGLE POINT ENERGY' in line:
                    data = line.strip().split()
                    qmenergy = float(data[-1]) #- calculate_nuclear_self_interaction(qmfile)
                    break
                else:
                    continue

        '''
        with open(str(qmfile + ".log")) as ifile:
            for line in ifile:
                match = re.search(r"^\s*LOEWDIN\s*ATOMIC\s*CHARGES", line, flags=re.MULTILINE) #match = re.search(r"^\s*CHELPG\s*Charges", line, flags=re.MULTILINE)  # or match = re.search(r"^\s*MULLIKEN\s*ATOMIC\s*CHARGES", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        break
                    for line in ifile:
                        match = re.search(
                            r"^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        match2 = re.search(
                            r"^\s*(\d+)\s+(\S+)\s+(\S+)", line, flags=re.MULTILINE
                        )
                        if match:
                            qm_corrdata.append(
                                [
                                    int(match.group(1)) + 1,
                                    match.group(2),
                                    float(match.group(4)),
                                ]
                            )
                        elif match2:
                            elem = str(match2.group(2)).rstrip(':')# if elem.endswith(':') else str(match2.group(2))
                            qm_corrdata.append(
                                [
                                    int(match2.group(1)) + 1,
                                    elem,
                                    float(match2.group(3)),
                                ]
                            )
                        else:
                            break
                    break
    if str(qmprog) == "SERENITY":
        if fde == 0:
            with open(str(qmfile + ".log")) as ifile:
                lines = ifile.readlines()
                for line in lines:
                    if 'Total Energy' in line:
                        data=line.strip().split()
                        qmenergy = float(data[-1])
                        break
                    else:
                        continue
            
                stop= False
                for line in range(len(lines)):
                    if 'Hirshfeld Population Analysis' in lines[line]:
                        for i in range(line+5, len(lines), 1):
                            if lines[i].strip() != "" and 'Time taken for task' not in lines[i+1]:
                                data = lines[i].strip().split()
                                if len(data) >= 4:
                                    qm_corrdata.append([
                                                int(data[0]),
                                                data[1],
                                                float(data[3]),
                                            ]
                                        )
                            else:
                                stop=True
                                break
                    if stop:
                        break
        elif fde == 1:
            with open(str(qmfile + ".log")) as ifile:
                lines = ifile.readlines()
                FSR_found = False
                hirshfeld_section_found = False
                FaT_energies_found = False

                for line in lines:
                    if tddft == 1 and 'Final Freeze-and-Thaw Energies' in line:
                        FaT_energies_found = True

                    elif 'Final SCF Results' in line and (not tddft == 1):
                        FSR_found = True

                    if FSR_found or FaT_energies_found:
                        if 'Total Supersystem Energy (active + all env.):' in line:
                            data = line.strip().split()
                            qmenergy = float(data[-1])  # Extract the energy value, which is the last element
                            break

                if FSR_found or FaT_energies_found:
                    for i in range(len(lines)):
                        if 'Hirshfeld Population Analysis' in lines[i]:
                            hirshfeld_section_found = True

                        if hirshfeld_section_found:
                            for j in range(i + 5, len(lines)):
                                if j + 1 < len(lines) and lines[j].strip() != "" and 'Time taken for task' not in lines[j + 1]:
                                    data = lines[j].strip().split()
                                    if len(data) >= 4:
                                        qm_corrdata.append([
                                            int(data[0]),  # Atom number
                                            data[1],       # Atom type
                                            float(data[3]) # Charge
                                        ])
                                else:
                                    hirshfeld_section_found = False
                                    break

                for i in range(len(qm_corrdata)):
                    qm_corrdata[i][0] = i + 1  

            '''with open(str(qmfile + ".log")) as ifile:
                lines = ifile.readlines()
                FSR_found = False
                for line in lines:
                    if 'Final SCF Results' in line:
                        FSR_found = True
                    if FSR_found:
                        if 'Total Supersystem Energy (active + all env.):' in line:
                            data=line.strip().split()
                            qmenergy = float(data[7])
                            break
                        else:
                            continue
                stop= False
                hirshfeld_section_found = False
                for line in range(len(lines)):
                    if 'Hirshfeld Population Analysis' in lines[line]:
                        hirshfeld_section_found = True
                    if hirshfeld_section_found:
                        for i in range(line+5, len(lines), 1):
                            if lines[i].strip() != "" and 'Time taken for task' not in lines[i+1]:
                                data = lines[i].strip().split()
                                if len(data) >= 4:
                                    qm_corrdata.append([
                                                int(data[0]),
                                                data[1],
                                                float(data[3]),
                                            ]
                                        )
                            elif lines[i].strip() == "" and 'Time taken for task' in lines[i+1] and 'Task No.' in lines[i+4]:
                                hirshfeld_section_found = False
                                stop=False
                                break
                            else:
                                hirshfeld_section_found = False
                                stop=True
                                break
                    if stop:
                        break
            for i in range(len(qm_corrdata)):
                qm_corrdata[i][0] = i+1'''

    if str(qmprog) == "TM":
        atoms = []
        with open(str(qmfile[:-4] + ".str")) as ifile:
            lines = ifile.readlines()
            for line in lines:
                data = line.strip().split()
                if len(data) == 4:
                    atoms.append(data[0])
                else:
                    continue

        with open(str(qmfile + ".log")) as ifile:
            lines = ifile.readlines()
            for line in lines:
                data = line.strip().split()
                if 'total energy' in line and (len(data) == 6):
                    qmenergy = float(data[-2])  # Extract the energy value, which is the last element
                    break
                else:
                    continue
            
            mk_section_found = False
            for i in range(len(lines)):
                if 'charges resulting from fit:' in lines[i]:
                    mk_section_found = True

                if mk_section_found:
                    count = 0
                    for j in range(i + 4, i + 4 + len(atoms)):
                        data = lines[j].strip().split()
                        count += 1
                        qm_corrdata.append([
                            count,                # Atom number
                            atoms[count-1],       # Atom type
                            float(data[-1])       # Charge
                        ])
                    mk_section_found = False
                else:
                    continue

    logger(logfile, "QM energy is " + str(float(qmenergy)) + " a.u..\n")
    return qmenergy, qm_corrdata

def get_mmenergy(edrname, qmmmInputs):
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    logfile = qmmmInputs.logfile

    mmenergy = 0.0
    logger(logfile, "Extracting MM energy.\n")
    p = subprocess.Popen(
        [
            prefix,
            "energy",
            "-f",
            edrname,
            "-o",
            str(edrname + ".xvg"),
            "-backup",
            "no",
        ],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    p.communicate(input=b"11\n\n")
    with open(str(edrname + ".xvg")) as ifile:
        for line in ifile:
            match = re.search(
                r"^    0.000000\s*([-]*\d+.\d+)\n", line, flags=re.MULTILINE
            )
            if match:
                mmenergy = float(match.group(1)) * 0.00038087988
                break
    logger(logfile, "MM energy is " + str(float(mmenergy)) + " a.u..\n")
    return mmenergy

def take_third(elem):
    import numpy as np
    out=np.sqrt(float(elem[2])*float(elem[2]))
    return out

def get_statechar(qmfile, qmmmInputs):
    qmprog = qmmmInputs.qmparams.program
    extra_string = qmmmInputs.qmparams.extra
    pcffile = qmmmInputs.pcffile
    logfile = qmmmInputs.logfile
    basedir = qmmmInputs.basedir
    logger(logfile, "Extracting state character.\n")
    char=[]
    if str(qmprog) == "G16":
        find=False
        textfile=open(str(qmfile + ".log"))
        ifile = textfile.readlines()
        nl=0
        for line in reversed(ifile):
            if not find:
                match = re.search(r"^\s*This\sstate\sfor\soptimization.*",line, flags=re.MULTILINE)
                if match:
                    find=True
            else:
                match=re.match(r"^\s*(\d+)\s\S{2}(\d+)\s*([-\d]\d+\.\d+).*",line, flags=re.MULTILINE)
                match2=re.match(r"^\s*(\d+)\s\S{2}(\d+)\s*(\d+\.\d+).*",line, flags=re.MULTILINE)
                if match:
                    char.append([match.group(1),match.group(2),match.group(3)])
                elif match2:
                    char.append([match2.group(1),match2.group(2),match2.group(3)])
                else:
                   break
        if char==[]:
             logger(logfile, "Couldn't find the character of the optimized excited state.\n")
             exit(1)
    textfile.close()
    char.sort(key=take_third,reverse=True)
    return char

def get_all_states(qmfile, qmmmInputs):
    extra_string = qmmmInputs.qmparams.extra
    logfile = qmmmInputs.logfile
    logger(logfile, "Extracting state character.\n")
    chars=[]
    if str(qmprog) == "G16":
        find=False
        textfile=open(str(qmfile + ".log"))
        ifile = textfile.readlines()
        nl=0
        for line in reversed(ifile):
            match = re.search(r"^\s*Excitation\senergies\sand\soscillator\sstrengths.*",line, flags=re.MULTILINE)
            if match:
                break
            else:
                nl+=1
        new_state=True
        for line in ifile[-nl-1:]:
            if new_state:
                char=[]
            match=re.match(r"^\s*(\d+)\s\S{2}(\d+)\s*([-\d]\d+\.\d+).*",line, flags=re.MULTILINE)
            match2=re.match(r"^\s*(\d+)\s\S{2}(\d+)\s*(\d+\.\d+).*",line, flags=re.MULTILINE)
            if match:
               char.append([match.group(1),match.group(2),match.group(3)])
               new_state=False
            elif match2:
               char.append([match2.group(1),match2.group(2),match2.group(3)])
               new_state=False
            elif not new_state:
                new_state=True
                char.sort(key=take_third,reverse=True)
                chars.append(char[0:3])
        if chars==[]:
            logger(logfile, "Couldn't find the character of the optimized excited state.\n")
            exit(1)
    return chars

def get_linkenergy_au(qm_corrdata, qmmmInputs):
    xyzq = qmmmInputs.xyzq
    linkcorrlist = qmmmInputs.linkcorrlist
    m1list = qmmmInputs.m1list
    m2list = qmmmInputs.m2list
    q1list = qmmmInputs.q1list
    qmmmtop = qmmmInputs.qmmmtop
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile
    pcffile = qmmmInputs.pcffile
    qmatomlist = qmmmInputs.qmatomlist

    linkenergy = 0.0
    m2charges = get_m2charges(xyzq, m1list, m2list)
    for element in linkcorrlist:
        z1 = 0.0
        v1 = []
        v2 = []
        if int(element[0]) in np.array(list(_flatten(m2list))).astype(int):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[0]):
                        z1 = float(m2charges[i][j])
                        v1 = [
                            xyzq[int(element[0]) - 1][0] / 0.52917721,
                            xyzq[int(element[0]) - 1][1] / 0.52917721,
                            xyzq[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
                if z1 != 0.0:
                    break
        elif int(element[0]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i][2])
                    v1 = [
                        xyzq[int(element[0]) - 1][0] / 0.52917721,
                        xyzq[int(element[0]) - 1][1] / 0.52917721,
                        xyzq[int(element[0]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[0]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v1 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z1 = float(xyzq[int(element[0]) - 1][3])
            v1 = [
                xyzq[int(element[0]) - 1][0] / 0.52917721,
                xyzq[int(element[0]) - 1][1] / 0.52917721,
                xyzq[int(element[0]) - 1][2] / 0.52917721,
            ]
        z2 = 0.0
        if int(element[1]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[1]):
                        z2 = float(m2charges[i][j])
                        v2 = [
                            xyzq[int(element[1]) - 1][0] / 0.52917721,
                            xyzq[int(element[1]) - 1][1] / 0.52917721,
                            xyzq[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
                if z2 != 0.0:
                    break
        elif int(element[1]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i][2])
                    v2 = [
                        xyzq[int(element[1]) - 1][0] / 0.52917721,
                        xyzq[int(element[1]) - 1][1] / 0.52917721,
                        xyzq[int(element[1]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[1]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v2 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z2 = float(xyzq[int(element[1]) - 1][3])
            v2 = [
                xyzq[int(element[1]) - 1][0] / 0.52917721,
                xyzq[int(element[1]) - 1][1] / 0.52917721,
                xyzq[int(element[1]) - 1][2] / 0.52917721,
            ]
        v12 = np.array(v1) - np.array(v2)
        dist = np.linalg.norm(v12)
        linkenergy += z1 * z2 / dist
    # now also all atoms in the corrdata list with the mod and linkcorr point charges
    # mod first. mod is charge in pcffile minus m2charge
    pcf = read_pcffile(pcffile)
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(m2list[i][j]) - 1][k]) / 0.52917721)
            curr_mod_charge = (
                float(float(pcf[int(m2list[i][j]) - 1][3])) - m2charges[i][j]
            )
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
    # now linkcorr. linkcorr are last m2*2 entries in pcf
    m2count = 0
    linkstart = len(pcf) - 2 * len(list(_flatten(m2list)))
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
            curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
            m2count += 1
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                linkenergy += z1 * curr_mod_charge / dist
    # now, add the correction of energy for the link atoms. currently only C-C bond cuts supported.
    for i in range(0, len(linkatoms)):
        v1 = [
            linkatoms[i][0] / 0.52917721,
            linkatoms[i][1] / 0.52917721,
            linkatoms[i][2] / 0.52917721,
        ]
        _flattened = list(_flatten(q1list))
        v2 = [
            xyzq[int(_flattened[i]) - 1][0] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][1] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][2] / 0.52917721,
        ]
        v12 = np.array(v2) - np.array(v1)
        dist = np.linalg.norm(v12)
    dist = dist * 0.7409471631
    energycorr = databasecorrection("ENERGY", "aminoacid_CACB", dist, qmmmInputs)
    linkenergy -= energycorr
    # sign inverted due to correction convention (subtracting)
    return linkenergy

def get_energy(qmfile, edrname, qmmmInputs):
    logfile = qmmmInputs.logfile
    qmprog = qmmmInputs.qmparams.program #Added by Nicola
    qmenergy, qm_corrdata = get_qmenergy(qmfile, qmmmInputs)
    mmenergy = get_mmenergy(str(edrname), qmmmInputs)
    if qmmmInputs.linkcorrlist:# and qmprog == "G16": #Added by Nicola
        linkcorrenergy = get_linkenergy_au(qm_corrdata, qmmmInputs)
    else:
        linkcorrenergy = 0.0
    basis = qmmmInputs.qmparams.basis
    #qmenergy -= linkcorrenergy
    methodstring = str(qmmmInputs.qmparams.method)
    if basis != "NONE":
        methodstring += str("/" + str(basis))
    logger(
        logfile,
        str(
            "Single point energy done. QM/MM energy is {:>20.10f} (QM, link atom corrected ".format(
                float(qmenergy-linkcorrenergy)
            )
            + methodstring
            + ") + {:>20.10f} (MM) = {:>20.10f} (a.u.)\n".format(
                float(mmenergy), float(qmenergy-linkcorrenergy) + float(mmenergy)
            )
        ),
    )

    return qmenergy, mmenergy, linkcorrenergy, qm_corrdata

#Forces
def transform_D_to_E(num):
    number = num.replace('D', 'E')
    number = float(number)
    return number

def get_qmforces_au(qmmmInputs):
    qmatomlist = qmmmInputs.qmatomlist
    m1list = qmmmInputs.m1list
    qmmmtop = qmmmInputs.qmmmtop
    qmprogram = qmmmInputs.qmparams.program
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    logfile = qmmmInputs.logfile
    if qmprogram == "SERENITY":
        fde = qmmmInputs.qmparams.fde

    qmforces = []
    qmonlyforcelist = []
    pcf_grad = []

    if qmprogram == "G16":
        insert = ""
        #2222
        if (int(curr_step) != 0): #and (qmmmInputs.qmmmparams.jobname == 'OPT'):  why should the jobname be relevant ??? #SIMON
            insert = str("." + str(curr_step))  #change curr_step-1 to curr_step   #SIMON WAR HIER AM WERK
        
        logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".gjf.log")+" and "+ str(jobname + insert + ".fort.7")+"\n")
        qmlogfile = str(jobname + insert + ".gjf.log")
        fortfile = str(jobname + insert + ".fort.7")

        with open(fortfile) as ifile:
            for line in ifile:
                match = re.search(
                    r"^\s*(\S*)\s*(\S*)\s*(\S*)", line, flags=re.MULTILINE
                )
                if match:
                    qmline = [
                        float(str(match.group(1)).replace("D", "e")) * -1.0,
                        float(str(match.group(2)).replace("D", "e")) * -1.0,
                        float(str(match.group(3)).replace("D", "e")) * -1.0,
                    ]
                    qmonlyforcelist.append(qmline)
        with open(qmlogfile) as i2file:
            for line in i2file:
                match = re.search(
                    r"^\s*Electrostatic\s*Properties\s*\(Atomic\s*Units\)",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    for line in i2file:
                        match = re.search(
                            r"^\s*\S+\s*[-]*\d+\.\d+\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)\s*([-]*\d+\.\d+)",
                            line,
                            flags=re.MULTILINE,
                        )
                        if match:
                            pcf_grad_line = [
                                float(match.group(1)),
                                float(match.group(2)),
                                float(match.group(3)),
                            ]
                            pcf_grad.append(pcf_grad_line)
                        match = re.search(r"^\s*Leave Link", line, flags=re.MULTILINE)
                        if match:
                            break
                    break
    if qmprogram == "ORCA":#Added by Nicola
        insert = ""
        #2222
        if (int(curr_step) != 0): #and (qmmmInputs.qmmmparams.jobname == 'OPT'):  why should the jobname be relevant ??? #SIMON
            insert = str("." + str(curr_step))  #change curr_step-1 to curr_step   #SIMON WAR HIER AM WERK
        
        logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".inp.log")+" and "+ str(jobname + insert + ".engrad")+" and "+ str(jobname + insert + ".pcgrad") + "\n")
        qmlogfile = str(jobname + insert + ".inp.log")
        engradfile = str(jobname + insert + ".engrad")
        pcgradfile = str(jobname + insert + ".pcgrad")
        with open(qmlogfile) as ifile:
            lines = ifile.readlines()
            cartesian_gradient_section_found = False
            for line in lines:
                if 'CARTESIAN GRADIENT' in line:
                    cartesian_gradient_section_found = True
                    continue
                
                if cartesian_gradient_section_found:
                    data = line.strip().split()
                    if len(data) == 6:
                        if data[2] == ':':
                            #if data[1] != 'Q':
                            qmonlyforcelist.append([float(data[3]) * -1.0, float(data[4]) * -1.0, float(data[5]) * -1.0]) 
                            #else:
                            #    pcf_grad.append([float(data[3]) * -1.0, float(data[4]) * -1.0, float(data[5]) * -1.0]) 
                        else:
                            continue
                    elif 'TIMINGS' in line:
                        cartesian_gradient_section_found = False
                        if qmonlyforcelist == []: #and pcf_grad == []:
                            logger(logfile,'No QM forces found, Exiting')
                            print('No QM forces found, Exiting')
                            exit(1)
                        #elif qmonlyforcelist == [] or pcf_grad == []:
                        #    logger(logfile,'Some QM forces may be missing')
                        break
                    else:
                        continue
        
        with open(pcgradfile) as i2file:
            lines = i2file.readlines()
            lines = lines[1:]
            for line in lines:
                data = line.strip().split()
                pcf_grad.append([float(data[0]) * -1.0, float(data[1]) * -1.0, float(data[2]) * -1.0])
    
    if qmprogram == "SERENITY":
        tddft = qmmmInputs.qmparams.tddft
        insert = ""
        #2222
        if (int(curr_step) != 0): #and (qmmmInputs.qmmmparams.jobname == 'OPT'):  why should the jobname be relevant ??? #SIMON
            insert = str("." + str(curr_step))  #change curr_step-1 to curr_step   #SIMON WAR HIER AM WERK

        logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".inp.log") + "\n")
        qmlogfile = str(jobname + insert + ".inp.log")
        with open(qmlogfile) as ifile:
            lines = ifile.readlines()
            cartesian_gradient_section_found = False
            pc_gradient_section_found = False
            #lines = ifile.readlines()
            #stop = False
            for line in lines:
                if 'Current Geometry Gradients (a.u.):' in line:
                    cartesian_gradient_section_found = True
                    pc_gradient_section_found = False
                    continue
                elif 'Active System:' in line:
                    cartesian_gradient_section_found = False
                    pc_gradient_section_found = False
                    continue

                elif 'Point Charge Gradients (a.u.):' in line:
                    cartesian_gradient_section_found = False
                    pc_gradient_section_found = True
                    if qmonlyforcelist == []: #and pcf_grad == []:
                        logger(logfile,'No QM forces found, Exiting')
                        print('No QM forces found, Exiting')
                        exit(1)
                elif 'Time taken for task' in line and pc_gradient_section_found:
                    if pcf_grad == []:
                        logger(logfile,'No PCF forces found')
                    break


                if cartesian_gradient_section_found:
                    data = line.strip().split()
                    if len(data) == 5:
                        qmonlyforcelist.append([float(data[2]) * -1.0, float(data[3]) * -1.0, float(data[4]) * -1.0]) 
                    else:
                        continue
                elif pc_gradient_section_found:
                    data = line.strip().split()
                    if len(data) == 5:
                        pcf_grad.append([float(data[2]) * -1.0, float(data[3]) * -1.0, float(data[4]) * -1.0]) 
                    else:
                        continue

                
                
                        '''

                    for i in range(line+2, len(lines), 1):
                        if lines[i].strip() != "" or 'Time taken for task' not in lines[i+1]:
                            data = lines[i].strip().split()
                            qmonlyforcelist.append([float(data[2]) * -1.0, float(data[3]) * -1.0, float(data[4]) * -1.0]) 
                        else:
                            stop = True
                            break
                if stop:
                    break
        
        logger(logfile, "Point charge gradients are still not implemented in serenity, so the qm_pcf_grad will be set to zero")
        with open('pointcharges.pc') as i2file:
            lines = i2file.readlines()
            num = int(lines[0])
            for i in range(num):
                pcf_grad.append(float(0))
            '''     

    if qmprogram == "TM":

        insert = ""
        #2222
        if (int(curr_step) != 0): #and (qmmmInputs.qmmmparams.jobname == 'OPT'):  why should the jobname be relevant ??? #SIMON
            insert = str("." + str(curr_step))  #change curr_step-1 to curr_step   #SIMON WAR HIER AM WERK
        
        logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".inp.log") + "\n")
        qmlogfile = str(jobname + insert + ".inp.log")
        atoms = []
        with open(str(jobname + insert + ".str")) as ifile:
            lines = ifile.readlines()
            for line in lines:
                data = line.strip().split()
                if len(data) == 4:
                    atoms.append(data[0])
                else:
                    continue

        with open(qmlogfile) as ifile:
            lines = ifile.readlines()
            cartesian_gradient_section_found = False
            pc_gradient_section_found = False
            for line in range(len(lines)):
                if '$grad          cartesian gradients' in lines[line]:
                    cartesian_gradient_section_found = True
                    pc_gradient_section_found = False

                elif '$point_charge_gradients' in lines[line]:
                    pc_gradient_section_found = True
                    cartesian_gradient_section_found = False

                if cartesian_gradient_section_found:
                    for i in range(line + len(atoms) + 2, len(lines)):
                        data = lines[i].strip().split()
                        qmonlyforcelist.append([float(data[0].replace('D', 'E')) * -1.0, float(data[1].replace('D', 'E')) * -1.0, float(data[2].replace('D', 'E')) * -1.0])
                        if '$end' in lines[i+1]:
                            cartesian_gradient_section_found = False
                            break
                        else:
                            continue

                elif pc_gradient_section_found:
                    for i in range(line + 1, len(lines)):
                        data = lines[i].strip().split()
                        if len(data) == 3:
                            pcf_grad.append([float(data[0].replace('D', 'E')) * -1.0, float(data[1].replace('D', 'E')) * -1.0, float(data[2].replace('D', 'E')) * -1.0])
                        else:
                            pc_gradient_section_found = False
                            break
                else:
                    continue
                        




    with open(qmmmtop) as ifile:
        for line in ifile:
            match = re.search(r"\[\s+moleculetype\s*\]", line)
            if match:
                for line in ifile:
                    match = re.search(r"\[\s+atoms\s*\]", line)
                    if match:
                        count = 0
                        qmcount = 0
                        m1count = 0
                        for line in ifile:
                            match = re.search(
                                r"^\s*(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-]*\d+\.*\d*)\s+(\d+\.*\d*)",
                                line,
                                flags=re.MULTILINE,
                            )
                            if (
                                match
                                and (
                                    int(match.group(1))
                                    not in np.array(qmatomlist).astype(int)
                                )
                                and (int(match.group(1)) not in np.array(m1list).astype(int))
                            ):
                                if qmprogram == "G16":
                                    curr_charge = float(match.group(7))
                                    qmforces.append(
                                        [
                                            pcf_grad[count][0] * curr_charge,
                                            pcf_grad[count][1] * curr_charge,
                                            pcf_grad[count][2] * curr_charge,
                                        ]
                                    )
                                    count += 1
                                elif qmprogram == "ORCA" or qmprogram == "SERENITY" or qmprogram == "TM":#Added by Nicola
                                    qmforces.append(pcf_grad[count])
                                    count += 1
                            elif match and int(match.group(1)) in np.array(
                                qmatomlist
                            ).astype(int):
                                qmforces.append(qmonlyforcelist[qmcount])
                                qmcount += 1
                            elif match and int(match.group(1)) in np.array(m1list).astype(
                                int
                            ):
                                qmforces.append(
                                    qmonlyforcelist[m1count + len(qmatomlist)]
                                )
                                m1count += 1
                            match = re.search(r"^\s*\n", line, flags=re.MULTILINE)
                            if match:
                                break
                        break
                break
    return qmforces

def get_mmforces_au(qmmmInputs):
    jobname = qmmmInputs.qmmmparams.jobname
    curr_step = qmmmInputs.qmmmparams.curr_step
    prefix =  qmmmInputs.pathparams.gmxpath + qmmmInputs.pathparams.gmxcmd
    logfile = qmmmInputs.logfile

    mmforces = []
    insert = ""
    if int(curr_step) != 0:
        insert = str("." + str(curr_step))
    logger(logfile,"Reading MM forces using file: "+str(jobname + insert + ".trr/.tpr/.tpr")+"\n")
    trrname = str(jobname + insert + ".trr")
    tprname = str(jobname + insert + ".tpr")
    xvgname = str(jobname + insert + ".xvg")
    p = subprocess.Popen(
        [
            prefix,
            "traj",
            "-fp",
            "-f",
            trrname,
            "-s",
            tprname,
            "-of",
            xvgname,
            "-xvg",
            "none",
            "-backup",
            "no",
        ],
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    p.communicate(input=b"0\n")
    with open(xvgname) as ifile:
        for line in ifile:
            forcelist = re.findall("\S+", line)
            count = 0
            mmforceline = []
            for i in range(1, len(forcelist)):
                mmforceline.append(float(forcelist[i]) * 2.0155295e-05)
                count += 1
                if count > 2:
                    count = 0
                    mmforces.append(mmforceline)
                    mmforceline = []
            break  # read only one line

    return mmforces

def get_linkforces_au(qm_corrdata,qmmmInputs):
    linkcorrlist = qmmmInputs.linkcorrlist
    qmatomlist = qmmmInputs.qmatomlist
    xyzq = qmmmInputs.xyzq
    pcffile = qmmmInputs.pcffile
    m1list = qmmmInputs.m1list
    m2list = qmmmInputs.m2list
    q1list = qmmmInputs.q1list
    linkatoms = qmmmInputs.linkatoms
    logfile = qmmmInputs.logfile

    linkforces = []
    # Force Coulomb: z1*z2*(distance along coord)/(distance between charges)**3
    for element in xyzq:  # this is just to count an entry for each atom!
        linkforces.append([0.0, 0.0, 0.0])
    m2charges = get_m2charges(xyzq, m1list, m2list)
    for element in linkcorrlist:
        z1 = 0.0
        v1 = []
        v2 = []
        if int(element[0]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[0]):
                        z1 = float(m2charges[i][j])
                        v1 = [
                            xyzq[int(element[0]) - 1][0] / 0.52917721,
                            xyzq[int(element[0]) - 1][1] / 0.52917721,
                            xyzq[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
                if z1 != 0.0:
                    break
        elif int(element[0]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i][2])
                    v1 = [
                        xyzq[int(element[0]) - 1][0] / 0.52917721,
                        xyzq[int(element[0]) - 1][1] / 0.52917721,
                        xyzq[int(element[0]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[0]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[0]):
                    z1 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v1 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z1 = float(xyzq[int(element[0]) - 1][3])
            v1 = [
                xyzq[int(element[0]) - 1][0] / 0.52917721,
                xyzq[int(element[0]) - 1][1] / 0.52917721,
                xyzq[int(element[0]) - 1][2] / 0.52917721,
            ]
        z2 = 0.0
        if int(element[1]) in _flatten(m2list):
            for i in range(0, len(m2list)):
                for j in range(0, len(m2list[i])):
                    if int(m2list[i][j]) == int(element[1]):
                        z2 = float(m2charges[i][j])
                        v2 = [
                            xyzq[int(element[1]) - 1][0] / 0.52917721,
                            xyzq[int(element[1]) - 1][1] / 0.52917721,
                            xyzq[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
                if z2 != 0.0:
                    break
        elif int(element[1]) in np.array(qmatomlist).astype(int):
            for i in range(0, len(qmatomlist)):
                if int(qmatomlist[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i][2])
                    v2 = [
                        xyzq[int(element[1]) - 1][0] / 0.52917721,
                        xyzq[int(element[1]) - 1][1] / 0.52917721,
                        xyzq[int(element[1]) - 1][2] / 0.52917721,
                    ]
                    break
        elif int(element[1]) in np.array(m1list).astype(int):
            for i in range(0, len(m1list)):
                if int(m1list[i]) == int(element[1]):
                    z2 = float(qm_corrdata[i + len(qmatomlist)][2])
                    v2 = [
                        linkatoms[i][0] / 0.52917721,
                        linkatoms[i][1] / 0.52917721,
                        linkatoms[i][2] / 0.52917721,
                    ]
                    break
        else:
            z2 = float(xyzq[int(element[1]) - 1][3])
            v2 = [
                xyzq[int(element[1]) - 1][0] / 0.52917721,
                xyzq[int(element[1]) - 1][1] / 0.52917721,
                xyzq[int(element[1]) - 1][2] / 0.52917721,
            ]
        v12 = np.array(v1) - np.array(v2)
        dist = np.linalg.norm(v12)

        for i in range(0, 3):
            linkforces[int(element[0]) - 1][i] += (
                z1 * z2 * v12[i] / (dist * dist * dist)
            )
            linkforces[int(element[1]) - 1][i] -= (
                z1 * z2 * v12[i] / (dist * dist * dist)
            )
    # now also all atoms in the corrdata list with the mod and linkcorr point charges
    # mod first. mod is charge in pcffile minus m2charge
    pcf = read_pcffile(pcffile)
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i])):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(m2list[i][j]) - 1][k]) / 0.52917721)
            curr_mod_charge = (
                float(float(pcf[int(m2list[i][j]) - 1][3])) - m2charges[i][j]
            )
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(qmatomlist[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
                    linkforces[int(m2list[i][j]) - 1][l] -= (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(m1list[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
                    linkforces[int(m2list[i][j]) - 1][l] -= (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
    m2count = 0
    linkstart = len(pcf) - 2 * len(list(_flatten(m2list)))
    for i in range(0, len(m2list)):
        for j in range(0, len(m2list[i]) * 2):
            curr_mod = []
            for k in range(0, 3):
                curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
            curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
            m2count += 1
            for k in range(0, len(qmatomlist)):
                v1 = [
                    xyzq[int(qmatomlist[k]) - 1][0] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][1] / 0.52917721,
                    xyzq[int(qmatomlist[k]) - 1][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(qmatomlist[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
            for k in range(0, len(linkatoms)):
                v1 = [
                    linkatoms[k][0] / 0.52917721,
                    linkatoms[k][1] / 0.52917721,
                    linkatoms[k][2] / 0.52917721,
                ]
                z1 = float(qm_corrdata[k + len(qmatomlist)][2])
                v12 = np.array(v1) - np.array(curr_mod)
                dist = np.linalg.norm(v12)
                for l in range(0, 3):
                    linkforces[int(m1list[k]) - 1][l] += (
                        z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                    )
    for i in range(0, len(linkatoms)):
        v1 = [
            linkatoms[i][0] / 0.52917721,
            linkatoms[i][1] / 0.52917721,
            linkatoms[i][2] / 0.52917721,
        ]
        _flattened = list(_flatten(q1list))
        v2 = [
            xyzq[int(_flattened[i]) - 1][0] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][1] / 0.52917721,
            xyzq[int(_flattened[i]) - 1][2] / 0.52917721,
        ]
        v12 = np.array(v2) - np.array(v1)
        dist = np.linalg.norm(v12) / 0.71290813568205
        u_v12 = rot.uvec(v12)
        dist = dist * 0.5282272551
        forcecorr = databasecorrection("FORCES", "aminoacid_CACB", dist, qmmmInputs)
        for j in range(0, 3):
            linkforces[int(_flattened[i]) - 1][j] += -u_v12[j] * forcecorr * 0.5
            linkforces[int(m1list[i]) - 1][j] += u_v12[j] * forcecorr * 0.5
    return linkforces

def read_forces(qm_corrdata,qmmmInputs):
    logfile = qmmmInputs.logfile
    active = qmmmInputs.active
    qmforces = []
    mmforces = []
    qmforces = get_qmforces_au(qmmmInputs)
    logger(logfile, str("QM forces read.\n"))
    mmforces = get_mmforces_au(qmmmInputs)
    logger(logfile, str("MM forces read.\n"))
    if qmmmInputs.linkcorrlist:
        linkcorrforces = get_linkforces_au(qm_corrdata, qmmmInputs)
    else:
        linkcorrforces = 0.0
    logger(logfile, str("Forces for link atom correction read.\n"))
    total_force = np.array(qmforces) + np.array(mmforces) - np.array(linkcorrforces)
    logger(logfile, str("Total forces obtained.\n"))
    if (qmmmInputs.qmmmparams.jobtype != "SINGLEPOINT") and (len(active) != 0):
        total_force=remove_inactive(total_force,active) #in SP case I don't have inactive atoms
        logger(logfile, str("Deleted forces of inactive atoms.\n"))
    return total_force

#update input
def make_opt_step(qmmmInputs):
    return 0

#scan inlcuded by SP (see nathaniel)
def scan_geo_dispvec(xyzq, stepsize, scan_atoms): #increment geometry
    dispvec = np.zeros((len(xyzq),3))
    coords = np.array(xyzq)[:, 0:3]
    #Bond length
    if len(scan_atoms) == 2 :
        a,b = scan_atoms
        coord_A, coord_B = (coords[int(a-1)], coords[int(b-1)])# array index -> -1
        line_vec = bmat.e_vector(coord_A, coord_B) * stepsize  # ba * stepsize 0.052917721 # 1 b,a.u. = 0.529177249 A 
        dispvec[int(b-1)] += line_vec

    #Angle
    elif len(scan_atoms) == 3 :
        a,b,c = scan_atoms
        coord_A, coord_B, coord_C = (coords[int(a-1)], coords[int(b-1)], coords[int(c-1)])# array index -> -1

        theta = stepsize * 180 / np.pi #degree -> radius
        #Dirty method
        bc_length = bmat.length(coord_B, coord_C)
        ba, bc = coord_B - coord_A, coord_B - coord_C
        b_normal = np.cross(ba, bc) # rotation axis
        bc_ = bmat.rot_mat(b_normal, theta) @ ba
        dispvec[int(c-1)] += bc_ - bc #cc_

    #Dihedral
    elif len(scan_atoms) == 4 :
        dihedral = bmat.angle(coords[scan_atoms[0]], coords[scan_atoms[1]], coords[scan_atoms[2]], coords[scan_atoms[3]])

    return dispvec

def scan_adj_force(forces, xyzq, scan_atoms, scan_flag='R') :
    coords = np.array(xyzq)[:, 0:3]

    if scan_flag == 'R': # bond  remove forces
        index = np.array([int(scan_atoms[0]-1), int(scan_atoms[1]-1)])# array index -> -1
        B_matrix = bmat.two_atom_B_matrix(coords[index[0]], coords[index[1]])
        fix_f = np.zeros(3) # init [0,0,0]

        for i in range(len(index)):
            fix_f += np.multiply(forces[index[i]], B_matrix[i])
        fix_f /= len(index)

        for i in range(len(index)):
            forces[index[i]] = forces[index[i]] - np.divide(fix_f,B_matrix[i])

    elif scan_flag == 'A': # angle
        index = np.array([int(scan_atoms[0]-1), int(scan_atoms[1]-1), int(scan_atoms[2]-1)])# array index -> -1
        B_matrix = bmat.three_atom_B_matrix(coords[index[0]], coords[index[1]], coords[index[2]])

    elif scan_flag == "D": # dihedral
        index = np.array([int(scan_atoms[0]-1), int(scan_atoms[1]-1), int(scan_atoms[2]-1), int(scan_atoms[3]-1)])# array index -> -1
        B_matrix = bmat.four_atom_B_matrix(coords[index[0]], coords[index[1]], coords[index[2]],coords[index[3]])


    return forces


#Job
def perform_sp(qmmmInputs):
    jobtype = qmmmInputs.qmmmparams.jobtype
    qmprogram = qmmmInputs.qmparams.program
    logfile = qmmmInputs.logfile
    curr_step = qmmmInputs.qmmmparams.curr_step

    logger(logfile, "Computing a single point.\n")
    logger(logfile, "Preparing QM and MM inputs:\n")
    if jobtype == "OPT":
        logger(logfile, "-----OPT %d---------------.\n"%curr_step)

    if qmprogram == "G16":
        #qm
        qmfile = make_g16_inp(qmmmInputs)
        logger(logfile, "Gaussian input file, %s, is ready.\n"%qmfile)
        print('-----Run g16file:%s---------------\n'%qmfile)
        run_g16(qmfile, qmmmInputs)

        #mm
        mmfile = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input file, %s, is ready.\n"%mmfile)      
        edrname = run_gmx(mmfile, qmmmInputs)
        
        logger(logfile, str("Reading energies.\n"))
        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(qmfile, edrname, qmmmInputs)
        total_energy = qmenergy + mmenergy - linkcorrenergy
        energies = (qmenergy, mmenergy, linkcorrenergy, total_energy)
        logger(logfile, str("Reading forces.\n"))
        total_force = read_forces(qm_corrdata,qmmmInputs)

        qmmmInputs.energies = energies
        qmmmInputs.forces = total_force

    elif qmprogram == "TM":
        #qm
        qmfile = make_TM_inp(qmmmInputs)
        logger(logfile, "Turbomole input file, %s, is ready.\n"%qmfile)
        print('-----Run Turbomole file:%s---------------\n'%qmfile)
        run_TM(qmfile, qmmmInputs)

        #mm
        mmfile = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input file, %s, is ready.\n"%mmfile)      
        edrname = run_gmx(mmfile, qmmmInputs)

        logger(logfile, str("Reading energies.\n"))
        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(qmfile, edrname, qmmmInputs)
        total_energy = qmenergy + mmenergy - linkcorrenergy
        energies = (qmenergy, mmenergy, linkcorrenergy, total_energy)

        logger(logfile, str("Reading forces.\n"))
        total_force = read_forces(qm_corrdata,qmmmInputs)

        qmmmInputs.energies = energies
        qmmmInputs.forces = total_force
        #logger(logfile,"Turbomole is not avalible currently.\n")
        #exit(0)

    elif qmprogram == "ORCA":#Added by Nicola
        #qm
        qmfile = make_orca_inp(qmmmInputs)
        logger(logfile, "Orca input file, %s, is ready.\n"%qmfile)
        print('-----Run orcafile:%s---------------\n'%qmfile)
        run_orca(qmfile, qmmmInputs)

        #mm
        mmfile = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input file, %s, is ready.\n"%mmfile)      
        edrname = run_gmx(mmfile, qmmmInputs)

        logger(logfile, str("Reading energies.\n"))
        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(qmfile, edrname, qmmmInputs)
        total_energy = qmenergy + mmenergy - linkcorrenergy
        energies = (qmenergy, mmenergy, linkcorrenergy, total_energy)

        logger(logfile, str("Reading forces.\n"))
        total_force = read_forces(qm_corrdata,qmmmInputs)

        qmmmInputs.energies = energies
        qmmmInputs.forces = total_force
        #logger(logfile,"ORCA is not avalible currently.\n")
        #exit(0)
    elif qmprogram == "SERENITY":
        #qm
        qmfile = make_serenity_inp(qmmmInputs)
        logger(logfile, "Serenity input file, %s, is ready.\n"%qmfile)
        print('-----Run Serenityfile:%s---------------\n'%qmfile)
        run_serenity(qmfile, qmmmInputs)

        #mm
        mmfile = make_gmx_inp(qmmmInputs)
        logger(logfile, "Gromacs input file, %s, is ready.\n"%mmfile)      
        edrname = run_gmx(mmfile, qmmmInputs)

        logger(logfile, str("Reading energies.\n"))
        qmenergy, mmenergy, linkcorrenergy, qm_corrdata = get_energy(qmfile, edrname, qmmmInputs)
        total_energy = qmenergy + mmenergy - linkcorrenergy
        energies = (qmenergy, mmenergy, linkcorrenergy, total_energy)

        logger(logfile, str("Reading forces.\n"))
        total_force = read_forces(qm_corrdata,qmmmInputs)

        qmmmInputs.energies = energies
        qmmmInputs.forces = total_force

    else:
        logger(logfile,"Unidentified QM program. Please check your -qm intput.\n")
        exit(0)

    # write a total force file
    if jobtype == "SCAN":
        logger(logfile, "Writing SCAN output.\n")

    elif jobtype == "SINGLEPOINT": # or curr_step == 0: #not nessesary ?
        logger(logfile, "Writing SP output.\n")
        write_output(energies, total_force, curr_step, energy_file="oenergy.txt", forces_file="oforces.txt")

    elif jobtype == "OPT":
        logger(logfile, "Writing OPT output.\n")

    elif jobtype == "NMA":
        logger(logfile, "Writing NMA output.\n")

 
def perform_opt(qmmmInputs):
    import copy
    #define done
    STEPLIMIT = 0
    FTHRESH = 1
    STEPSIZE = 2 

    logfile = qmmmInputs.logfile
    propagator = qmmmInputs.qmmmparams.propagater
    curr_step = qmmmInputs.qmmmparams.curr_step
    maxcycle = qmmmInputs.qmmmparams.maxcycle
    f_thresh = qmmmInputs.qmmmparams.f_thresh
    optlastonly = qmmmInputs.qmmmparams.optlastonly
    stepsize = qmmmInputs.qmmmparams.initstep
    jobname = qmmmInputs.qmmmparams.jobname
    jobtype = qmmmInputs.qmmmparams.jobtype
    qmprog = qmmmInputs.qmparams.program

    gro = qmmmInputs.gro
    xyzq = qmmmInputs.xyzq
    new_xyzq = []
    count = 0 #for each calculation counting for archiving
    count_trash=0 #just for SImon

    if curr_step == 0:
        #First calculation 
        perform_sp(qmmmInputs)
        old_qmmmInputs = copy.deepcopy(qmmmInputs)
        if jobtype == "SCAN" :
              write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step, energy_file="oenergy_%s.txt"%jobname, forces_file="oforces_%s.txt"%jobname)
        else:                                                                          
              write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step)

        #write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step)
        #Store 1st result
        curr_energy = qmmmInputs.energies[-1] #total energy
        total_force = qmmmInputs.forces
    else:
        curr_energy = read_oenergy(curr_step) #total energy
        total_force = read_oforces(curr_step)
        #Store previous result
        old_qmmmInputs = copy.deepcopy(qmmmInputs)
    
    last_forces = []
    clean_force = make_clean_force(total_force)
    maxforce = 0.0
    done = STEPLIMIT
    improved = True

    for element in _flatten(clean_force):
        if abs(float(element)) > abs(maxforce):
            maxforce = float(element)
    if abs(maxforce) < float(f_thresh):
        logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
        done = FTHRESH    
    else:
        logger(logfile, "Max force not below threshold. Continuing.\n")
        archive_first_step(jobname)                 #Simon

        # init BFGS hessian: initial guess for identity hessian
        if propagator == "BFGS":
            init_hessian = np.eye(3 * len(xyzq))
            np.savetxt("bfgs_hessian.txt", init_hessian)

        ############## start optimization loop ##############
        while not done and count <= maxcycle:
            #Store previous result
            #old_qmmmInputs = copy.deepcopy(qmmmInputs)

            
            #Add step first (if not improve, reject -1)
            qmmmInputs.qmmmparams.curr_step += 1
            curr_step = qmmmInputs.qmmmparams.curr_step
            
	    

            #Prepare new input            
            new_gro = str(jobname + "." + str(curr_step) + ".g96")
            new_pcffile = str(jobname + "." + str(curr_step) + ".pointcharges")
            #qmmmInputs.gro = new_gro
            #qmmmInputs.pcffile = new_pcffile
            logger(logfile,"using total force with max force: "+str(maxforce)+" for the displacement\n")
            if qmmmInputs.qmmmparams.jobtype == "SCAN" :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile, True, qmmmInputs.scan_atoms)
            else :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile)
            #    write_dispvec(dispvec, curr_step, count_trash) Simon implemented this to check if the dispvec is correct!  
            make_g96_inp(dispvec, gro, new_gro, logfile)
            qmmmInputs.gro = new_gro
            qmmmInputs.pcffile = new_pcffile
            qmmm_prep(qmmmInputs)
            
            #Run SP
            perform_sp(qmmmInputs) #make g16 & gmx input

            #Store new result
            curr_energy = qmmmInputs.energies[-1] #total energy
            last_energy = old_qmmmInputs.energies[-1]
            total_force = qmmmInputs.forces
            last_forces = old_qmmmInputs.forces

            ############## start energy check & update stepsize ############## 
            if curr_energy > last_energy:
                logger(logfile, "Rejected one optimization step due to energy increasing. Trying again, with smaller step.\n")
                count_trash+=1

                improved = False
                stepsize *= 0.2
                logger(logfile, "Reduced stepsize = "+str(stepsize)+"\n")
                #Simon wants to archive the trash stuff, please delete this if you are not Simon
                archive_trash(jobname, curr_step, count_trash)

                #remove files                
                insert      = str("." + str(curr_step))
                trrname     = str(jobname + insert + ".trr")
                tprname     = str(jobname + insert + ".tpr")
                gmxlogname  = str(jobname + insert + ".gmx.log")
                edrname     = str(jobname + insert + ".edr")
                xvgname     = str(jobname + insert + ".edr.xvg")
                g16name     = str(jobname + insert + ".gjf.log")
                fortname    = str(jobname + insert + ".fort.7")
                pcfname     = str(jobname + insert + ".pointcharges")
                gjfname     = str(jobname + insert + ".gjf")
                inpname     = str(jobname + insert + ".inp")
                outname     = str(jobname + insert + ".inp.log")
                engradname  = str(jobname + insert + ".engrad")
                xyzname     = str(jobname + insert + "*" + ".xyz")
                basisname   = str(jobname + insert + ".basis")
                auxbasisname= str(jobname + insert + ".auxbasis")
                strname     = str(jobname + insert + ".str")
                coordname   = str(jobname + insert + ".coord")
                #mosname     = str(jobname + insert + ".mos")
                subprocess.call("rm %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s"%(trrname,tprname,gmxlogname,edrname,xvgname,g16name,fortname,pcfname,gjfname,inpname,outname,engradname,xyzname,basisname,auxbasisname,coordname,strname), shell=True)

            else:
                count_trash=0
                stepsize *= 1.2
                improved = True

                #archive remove previous
                logger(logfile, "Archive previous files\n")
                archive(jobname, curr_step-1) #archive previous


            ############## end energy check & update stepsize ##############
            
            ############## start force check ##############
            clean_force = make_clean_force(total_force)
            maxforce = 0.0
            
            for element in _flatten(clean_force):
                if abs(float(element)) > abs(maxforce):
                    maxforce = float(element)
            if abs(maxforce) < float(f_thresh):
                logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
                done = FTHRESH
                break
            ############## end force check ##############

            ############## stepsize check ##############

            if float(stepsize) < 1e-6: #0.000001 a.u.
                done = STEPSIZE
                logger(
                    logfile,
                        ("Step became lower than 0.000001 a.u., optimization is considered done for now. " + 
                         "This is the best we can do unless reaching unacceptable numerical noise levels.\n"),
                )
                break
            ############## end stepsize check ##############

            if improved:
                #Update for BFGS              
                xyzq = old_qmmmInputs.xyzq
                new_xyzq = qmmmInputs.xyzq
                #Store previous 
                old_qmmmInputs = copy.deepcopy(qmmmInputs)
                count += 1
                if jobtype == "SCAN" :
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, energy_file="oenergy_%s.txt"%jobname, forces_file="oforces_%s.txt"%jobname)
                else:
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step)
                gro = qmmmInputs.gro            #SIMON
                logger(logfile, "Due to the decrease of the energy, the structure "+str(gro)+" will be used from now on.\n")
            else:
                #rejected and use previous
                qmmmInputs.qmmmparams.curr_step -= 1
                qmmmInputs = copy.deepcopy(old_qmmmInputs)
                total_force = qmmmInputs.forces
                clean_force = make_clean_force(total_force)
                maxforce = 0.0
                for element in _flatten(clean_force):
                    if abs(float(element)) > abs(maxforce):
                        maxforce = float(element)
                logger(logfile,"Due to the rejection use maximum force: "+str(maxforce)+"\n")
                logger(logfile,"Due to the increase of the energy, the used structure remains "+str(gro)+".\n")


        ############## end optimization loop ##############
    #Remain first/last result 
    '''
    if qmmmInputs.qmmmparams.optlastonly == "YES":
        logger(logfile, "Remain last result. Remove other results.\n")
        #remove first
        remove_files(jobname, 0, True)
        for steps in range(count-1):
            filename = jobname+('.%d*'%(steps+1))
            subprocess.call("rm %s"%(filename), shell=True)
    '''
    # Remove the last but failed with criteria files         
    subprocess.call("rm %s"%str(jobname + '.' + str(count+1) + '*'), shell=True)
    
    # opt status 
    if done == STEPLIMIT:
        logger(logfile, "Optimization canceled due to step limit.\n")
        curr_step = qmmmInputs.qmmmparams.curr_step
    elif done == FTHRESH:
        logger(logfile, "Optimization finished due to force threshold.\n")
        curr_step = qmmmInputs.qmmmparams.curr_step
    elif done == STEPSIZE:
        logger(logfile, "Optimization finished due to step size.\n")
        curr_step = qmmmInputs.qmmmparams.curr_step
        qmmmInputs.qmmmparams.curr_step -= 1
        curr_step = qmmmInputs.qmmmparams.curr_step
    old_g96 = str(jobname + "." + str(curr_step) + ".g96")                                      #SP
    #final_gro = "opt.opt.gro"                                                                    
    #g96_to_gro(old_g96,final_gro,logfile)     # gro corrupted...... if no opt needed old_g96 is wrong .... if curr_step=0 should be without ".curr_step"..... 
    #logger(logfile,"Final geometry written to opt.opt.gro.")


def perform_opt_root(qmmmInputs):
    import copy
    #define done
    STEPLIMIT = 0
    FTHRESH = 1
    STEPSIZE = 2 

    logfile = qmmmInputs.logfile
    propagator = qmmmInputs.qmmmparams.propagater
    curr_step = qmmmInputs.qmmmparams.curr_step
    maxcycle = qmmmInputs.qmmmparams.maxcycle
    f_thresh = qmmmInputs.qmmmparams.f_thresh
    optlastonly = qmmmInputs.qmmmparams.optlastonly
    stepsize = qmmmInputs.qmmmparams.initstep
    jobname = qmmmInputs.qmmmparams.jobname
    qmprog = qmmmInputs.qmparams.program
    root = qmmmInputs.qmparams.extra.split("root=")[1].split(",")[0]

    gro = qmmmInputs.gro
    xyzq = qmmmInputs.xyzq
    new_xyzq = []
    count = 0 #for each calculation counting for archiving

    if curr_step == 0:
        #First calculation 
        perform_sp(qmmmInputs)
        old_qmmmInputs = copy.deepcopy(qmmmInputs)
        write_output(qmmmInputs.energies, qmmmInputs.forces, curr_step)
        qmmmInputs.statechar=get_statechar(qmfile, qmmmInputs)
        write_statechar(qmmmInputs.statechar, curr_step)
        #Store 1st result
        curr_energy = qmmmInputs.energies[-1] #total energy
        total_force = qmmmInputs.forces
    else:
        curr_energy = read_oenergy(curr_step) #total energy
        total_force = read_oforces(curr_step)

    last_forces = []
    clean_force = make_clean_force(total_force)
    maxforce = 0.0
    done = STEPLIMIT
    improved = True

    for element in _flatten(clean_force):
        if abs(float(element)) > abs(maxforce):
            maxforce = float(element)
    if abs(maxforce) < float(f_thresh):
        logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
        done = FTHRESH    
    else:
        logger(logfile, "Max force not below threshold. Continuing.\n")

        # init BFGS hessian: initial guess for identity hessian
        if propagator == "BFGS":
            init_hessian = np.eye(3 * len(xyzq))
            np.savetxt("bfgs_hessian.txt", init_hessian)

        ############## start optimization loop ##############
        while not done and count <= maxcycle:
            #Store previous result
            old_qmmmInputs = copy.deepcopy(qmmmInputs)

            
            #Add step first (if not improve, reject -1)
            qmmmInputs.qmmmparams.curr_step += 1
            curr_step = qmmmInputs.qmmmparams.curr_step
            
            #Prepare new input            
            new_gro = str(jobname + "." + str(curr_step) + ".g96")
            new_pcffile = str(jobname + "." + str(curr_step) + ".pointcharges")
            qmmmInputs.gro = new_gro
            qmmmInputs.pcffile = new_pcffile
            if qmmmInputs.qmmmparams.jobtype == "SCAN" :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile, True)
            else :
                dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, curr_step, logfile)
            make_g96_inp(dispvec, gro, new_gro, logfile)
            qmmm_prep(qmmmInputs)

            #Run SP
            perform_sp(qmmmInputs) #make g16 & gmx input

            #Store new result
            curr_energy = qmmmInputs.energies[-1] #total energy
            last_energy = old_qmmmInputs.energies[-1]
            total_force = qmmmInputs.forces
            last_forces = old_qmmmInputs.forces
            qmmmInputs.statechar=get_statechar(qmfile, qmmmInputs)
            curr_statechar=qmmmInputs.statechar
            last_statechar=old_qmmmInputs.statechar

            ############## start energy check & update stepsize ############## 
            if curr_energy > last_energy:
                logger(logfile, "Rejected one optimization step due to energy increasing. Trying again, with smaller step.\n")
                improved = False
                stepsize *= 0.2

                #remove files                
                insert = str("." + str(curr_step))
                trrname = str(jobname + insert + ".trr")
                tprname = str(jobname + insert + ".tpr")
                gmxlogname = str(jobname + insert + ".gmx.log")
                edrname = str(jobname + insert + ".edr")
                xvgname = str(jobname + insert + ".edr.xvg")
                g16name = str(jobname + insert + ".gjf.log")
                fortname = str(jobname + insert + ".fort.7")
                subprocess.call("rm %s %s %s %s %s %s %s"%(trrname,tprname,gmxlogname,edrname,xvgname,g16name,fortname), shell=True)
                
            elif (curr_statechar[0][0:2] == last_statechar[0][0:2] or curr_statechar[0][0:2] == last_statechar[1][0:2]) and (curr_statechar[1][0:2] == last_statechar[0][0:2] or curr_statechar[1][0:2] == last_statechar[1][0:2]):
                stepsize *= 1.2
                improved = True

                #archive remove previous
                logger(logfile, "Archive previous files\n")
                archive(jobname, curr_step-1) #archive previous
            else:
                logger(logfile, "Identified a root flipping. Changing the root\n")
                statechars=get_all_states(qmfile, qmmmInputs)
                for nitem,item in enumerate(statechars):
                    item_reduced=[]
                    for element in item:
                        item_reduced.append([element[0],element[1]])
                    if (last_statechar[0][0:2] in item_reduced) and (last_statechar[1][0:2] in item_reduced):
                        logger(logfile, "The state that you are looking for is state "+str(nitem)+"\n")
                        root = qmmmInputs.qmparams.extra.split("root=")[1].split(",")[0]
                        logger(logfile, str(root)+" --> "+str(nitem)+"\n")
                        qmmmInputs.qmparams.extra=qmmmInputs.qmparams.extra.split("root=")[0]+"root="+str(nitem)+","+qmmmInputs.qmparams.extra.split("root=")[1].split(",")[1]
                        break
            ############## end energy check & update stepsize ##############
            
            ############## start force check ##############
            clean_force = make_clean_force(total_force)
            maxforce = 0.0
            
            for element in _flatten(clean_force):
                if abs(float(element)) > abs(maxforce):
                    maxforce = float(element)
            if abs(maxforce) < float(f_thresh):
                logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
                done = FTHRESH
                break
            ############## end force check ##############

            ############## stepsize check ##############

            if float(stepsize) < 1e-6: #0.000001 a.u.
                done = STEPSIZE
                logger(
                    logfile,
                        ("Step became lower than 0.000001 a.u., optimization is considered done for now. " + 
                         "This is the best we can do unless reaching unacceptable numerical noise levels.\n"),
                )
                break
            ############## end stepsize check ##############

            if improved:
                #Update for BFGS              
                xyzq = old_qmmmInputs.xyzq
                new_xyzq = qmmmInputs.xyzq
                #Store previous 
                old_qmmmInputs = copy.deepcopy(qmmmInputs)
                count += 1
                if jobtype == "SCAN":
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, energy_file="oenergy_%s.txt"%jobname, forces_file="oforces_%s.txt"%jobname)
                else:
                    write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step)
            else:
                #rejected and use prvious
                qmmmInputs.qmmmparams.curr_step -= 1
                qmmmInputs = copy.deepcopy(old_qmmmInputs)


        ############## end optimization loop ##############
    #Remain first/last result 
    if qmmmInputs.qmmmparams.optlastonly == "YES":
        logger(logfile, "Remain last result. Remove other results.\n")
        #remove first
        remove_files(jobname, 0, True)
        for steps in range(count-1):
            filename = jobname+('.%d*'%(steps+1))
            subprocess.call("rm %s"%(filename), shell=True)
    else:
        for steps in range(count):
            remove_files(jobname,steps)

    # Remove the last but failed with criteria files         
    subprocess.call("rm %s"%str(jobname + '.' + str(count+1) + '*'), shell=True)
    
    # opt status 
    if done == STEPLIMIT:
        logger(logfile, "Optimization canceled due to step limit.\n")
    elif done == FTHRESH:
        logger(logfile, "Optimization finished due to force threshold.\n")
    elif done == STEPSIZE:
        logger(logfile, "Optimization finished due to step size.\n")
    last_step=curr_step-1
    old_g96 = str(jobname + "." + str(last_step) + ".g96")
    #final_gro = "opt.opt.gro"
    #g96_to_gro(old_g96,final_gro,logfile) #see comment above SP
    #logger(logfile,"Final geometry written to opt.opt.gro.")

def perform_scan(qmmmInputs):
    logfile = qmmmInputs.logfile
    #qmmmInputs.qmmmparams.propagator = "BFGS"
    r_array, a_array, d_array = qmmmInputs.qmmmparams.scan
    
    #Check scan input
    if (len(r_array)+len(a_array)+len(d_array)) == 0:
        logger(logfile, "No scan coordinates are found. Please check the scan file input.\n")
    else:
        logger(logfile, 'Read %d scan coordinates.\n'%(len(r_array)+len(a_array)+len(d_array)))
        origin_qmmmInputs = copy.deepcopy(qmmmInputs)
    
    #Check point
    check_flag = False
    scan_step = qmmmInputs.qmmmparams.scan_step
    if scan_step > 0 :
        check_flag = True
        logger(logfile, "Continue scan step%d\n"%qmmmInputs.qmmmparams.scan_step)

    #Single scan
    #Bond length scan
    if len(r_array) > 0:
        r_shape = r_array.shape
        if not check_flag :
            subprocess.call("mkdir scanR", shell=True)
        
        # Bond length scan loop
        for i in range(r_shape[0]):
            scan_atoms = (r_array[i][0], r_array[i][1])
            stepsize = r_array[0][-2]
            steps = r_array[0][-1]
            logger(logfile, "Scanning bond length between atom%d-%d\n"%scan_atoms)
            logger(logfile, "Scanned stepsize: %.3f, scan %d steps\n"%(stepsize, steps))
            
            origin_gro = origin_qmmmInputs.gro

            direc = "scanR/R%d-%d"%scan_atoms
            if not check_flag : 
                subprocess.call("mkdir %s"%direc, shell=True)
            
            enegy_arr = []
            # Bond length step optimization loop
            for j in range(int(steps)):
                #check point
                if (scan_step > (j+1)) and (scan_step != 0) :
                    continue

                subdirec = direc + "/step%d"%(j+1)
                if not check_flag :
                    subprocess.call("mkdir %s"%subdirec, shell=True)

                qmmmInputs = copy.deepcopy(origin_qmmmInputs)

            #    logger(logfile, "Start bond length scan step%d...\n"%(j+1))
                logger(logfile, "------ scan step%d ------\n"%(j+1))            # for clarity in logfile AJ

                logger(logfile, "Create new scanned geometry with increament %.3f\n"%(stepsize*(j+1)) )
                dispvec = scan_geo_dispvec(qmmmInputs.xyzq, (stepsize*(j+1)), scan_atoms)
                new_gro = "scanR%d-%d"%scan_atoms + "_%d.g96"%(j+1)
                make_g96_inp(dispvec, origin_gro, new_gro, logfile)
                
                qmmmInputs.gro = new_gro
                qmmmInputs.qmmmparams.jobname = "scanR%d-%d"%scan_atoms + "_%d"%(j+1)
                qmmmInputs.scan_atoms = scan_atoms
                
                if not check_flag :
                    subprocess.call("mkdir %s"%subdirec, shell=True)
                    qmmmInputs.qmmmparams.curr_step = 0

                else:
                    check_flag = False

                logger(logfile, ("Run optimization of scanR%d-%d"%scan_atoms + "_%d\n"%(j+1)) )
                perform_opt(qmmmInputs)

                filename = "scanR%d-%d"%scan_atoms + "_%d"%(j+1)
                subprocess.call("mv %s* %s"%(filename, subdirec), shell=True) 

                #read last energies from each optimization
                step_energies = read_oenergy(-1, filename = 'oenergy_%s.txt'%filename)
                step_forces = read_oforces(-1, filename = 'oforces_%s.txt'%filename)
                write_scan(step_energies, step_forces, (j+1))



                logger(logfile, "End bond length scan step%d...\n"%(j+1))   
    
            if qmmmInputs.qmmmparams.optlastonly == "Minimal":
                subprocess.call("rm oenergy_*", shell=True)
                subprocess.call("rm oforces_*", shell=True) 

    """
    #Angle scan
    if len(a_array) > 0:
        a_shape = a_array.shape
        print("a_shape", a_shape, '\n')
        subprocess.call("mkdir scanA", shell=True)
        for i in range(a_shape[0]):
            scan_atoms = (a_array[i][0], a_array[i][1], a_array[i][2])
            stepsize = a_array[0][-2]
            steps = a_array[0][-1]

            logger(logfile, "Scanning angle between atom%d-%d-%d\n"%scan_atoms)
            logger(logfile, "Scanned stepsize: %.3f degree, scan %d steps\n"%(stepsize, steps))
            #subprocess.call("mkdir scanA/A%d-%d-%d"%scan_atoms, shell=True)

            origin_gro = origin_qmmmInputs.gro

            # Angle step loop
            for j in range(int(steps)):
                qmmmInputs = copy.deepcopy(origin_qmmmInputs)
                logger(logfile, "Start angle scan step%d...\n"%(j+1))

                logger(logfile, "Create new scanned geometry with increament %.3f degree\n"%(stepsize*(j+1)) )
                dispvec = scan_geo_dispvec(qmmmInputs.xyzq, (stepsize*(j+1)), scan_atoms)
                new_gro = "scanA%d-%d-%d"%scan_atoms + "_%d.g96"%(j+1)
                make_g96_inp(dispvec, origin_gro, new_gro, logfile)

                qmmmInputs.gro = new_gro
                qmmmInputs.qmmmparams.jobname = "scanA%d-%d-%d"%scan_atoms + "_%d"%(j+1)
                qmmmInputs.qmmmparams.curr_step = 0
                qmmmInputs.scan_atoms = scan_atoms

                logger(logfile, ("Run optimization of scanA%d-%d-%d"%scan_atoms + "_%d\n"%(j+1)) )
                #perform_opt(qmmmInputs)

                logger(logfile, "End angle scan step%d...\n"%(j+1))

    #Dihedral angle scan
    if len(d_array) > 0:
        d_shape = d_array.shape
        print("d_shape", d_shape, '\n')
        #subprocess.call("mkdir scanD", shell=True)
        for i in range(d_shape[0]):
            scan_atoms = (d_array[i][0], d_array[i][1], d_array[i][2], d_array[i][3])
            stepsize = d_array[0][-2]
            steps = d_array[0][-1]

            logger(logfile, "Scanning dihedral angle between atom%d-%d-%d-%d\n"%scan_atoms)
            logger(logfile, "Scanned stepsize: %f, scan %d steps\n"%(stepsize, steps))
            #subprocess.call("mkdir scanD/D%d-%d-%d-%d"%scan_atoms, shell=True)
    """

    #Original: 0 Step
    #origin_qmmmInputs.qmmmparams.jobtype = "OPT"
    #perform_opt(origin_qmmmInputs)



def perform_nma(qmmmInputs):
    logfile = qmmmInputs.qmmmparams.logfile
    basedir = qmmmparams.basedir
    active = qmmmInputs.active
    jobname = qmmmInputs.qmmmparams.jobname

    logger(logfile, "------This will be a numerical) normal mode analysis.------\n")
    logger(
        logfile,
        "Generating a numerical Hessian for the active region using a displacement step of %f a.u.\n"%qmmmInputs.qmmmparams.disp,
    )

    logger(logfile, "Will require %d single point calculations!\n"%(len(active)*6+1))

    perform_sp(qmmmInputs)

    start_energy = qmmmInputs.energies[-1]
    start_forces = qmmmInputs.forces
    start_grad = np.array(start_forces) * -1.0

    hessian_xyz_full = []

    for curr_atom in active:
        grad_deriv_vec = nma_stuff.get_xyz_2nd_deriv(qmmmInputs, curr_atom, start_energy, start_forces)
        hessian_xyz_full.extend(grad_deriv_vec)
    prep_hess = nma_stuff.prepare_hess(hessian_xyz_full, active)

    for i in range(0, len(prep_hess[0]) - 6):
        evals.extend([float(1000.0)])
    
    nma_stuff.write_pseudofchk_file(
        jobname, evals, prep_hess, prep_hess, active, qmmmtop, logfile, xyzq
    )  # using prep_hess as pseudo-nm_matrix since we do not know yet
    
    logger(logfile, "Wrote pseudofchk (G03 format).\n")
    
    write_hess.hes_xyz_fchk(
        str(jobname + ".pseudofchk"), str(jobname + ".hess")
    )
    
    logger(logfile, "Wrote orca format .hess file.\n")
    evals, nm_matrix = nma.nma_3Nminus6dof_asfunction(
        str(jobname + ".hess"), basedir
    )
    print(nma_stuff.log_nma(
        qmmminfo, logfile, evals, nm_matrix, active, qmmmtop, xyzq, prep_hess
    ))

    return 0


if __name__ == "__main__":
    print("This file serves as a library for gmx2qmmm-related functions.")
    print("Do not execute directly. Use gmx2qmmm instead.")
