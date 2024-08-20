#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Run Singlepoint Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-12'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import math
import sqlite3
import sys
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators._helper import filter_xyzq, _flatten

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Singlepoint():

    '''
    This Class Performs A Singlepoint Calculation
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, basedir) -> None:
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        input_dict: dict -> Dictionary Of The User Parameter Input \\
        class_system: class -> Class Object Of The System\\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''
        # XX AJ we need to make a decision where we want to distinguish according to the QM software. Do we want all QM preparations in one file or one file per QM software?

        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.basedir = basedir

        #   XX AJ check how to deal with nma flag later
        self.nmaflag = 0
        
        # prepare QM input depending on software
        self.str_inputfile_qm = self.make_qm_input()

        # run QM calculation
        self.run_g16()

        # prepare MM input
        self.make_gmx_inp()
        # run MM calculation
        self.run_gmx()

        #   Read Energies
        # XX AJ rename variables
        self.qmenergy, self.mmenergy, self.linkcorrenergy, self.qm_corrdata = self.get_energy()
        self.total_energy = self.qmenergy + self.mmenergy - self.linkcorrenergy
        self.energies = (self.qmenergy, self.mmenergy, self.linkcorrenergy, self.total_energy)

        #   Read Forces
        self.total_force = self.read_forces()

        #   Write Output File
        if self.dict_input_userparameters['jobtype'] == 'SINGLEPOINT':
            self.write_output(energy_file="oenergy.txt", forces_file="oforces.txt")
        

    def get_energy(self):
        qmenergy, qm_corrdata = self.get_qmenergy()
        mmenergy = self.get_mmenergy()
        if self.system.linkcorrlist:
            linkcorrenergy = self.get_linkenergy_au(qm_corrdata)
        else:
            linkcorrenergy = 0.0
        methodstring = str(self.dict_input_userparameters['method'])
        if self.dict_input_userparameters['basis'] != "NONE":
            methodstring += str("/" + str(self.dict_input_userparameters['basis']))

        return qmenergy, mmenergy, linkcorrenergy, qm_corrdata
    
    
    def get_linkenergy_au(self, qm_corrdata):

        linkenergy = 0.0
        m2charges = self.get_m2charges()
        for element in self.system.linkcorrlist:
            z1 = 0.0
            v1 = []
            v2 = []
            if int(element[0]) in np.array(list(_flatten(self.system.list_atoms_m2))).astype(int):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[0]):
                            z1 = float(m2charges[i][j])
                            v1 = [
                                self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z1 != 0.0:
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[0]):
                        z1 = float(qm_corrdata[i][2])
                        v1 = [
                            self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[0]):
                        z1 = float(qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v1 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z1 = float(self.system.array_xyzq_current[int(element[0]) - 1][3])
                v1 = [
                    self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                ]
            z2 = 0.0
            if int(element[1]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[1]):
                            z2 = float(m2charges[i][j])
                            v2 = [
                                self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z2 != 0.0:
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[1]):
                        z2 = float(qm_corrdata[i][2])
                        v2 = [
                            self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[1]):
                        z2 = float(qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v2 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z2 = float(self.system.array_xyzq_current[int(element[1]) - 1][3])
                v2 = [
                    self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                ]
            v12 = np.array(v1) - np.array(v2)
            distance = np.linalg.norm(v12)
            linkenergy += z1 * z2 / distance
        # now also all atoms in the corrdata list with the mod and linkcorr point charges
        # mod first. mod is charge in pcffile minus m2charge
        # XX AJ the read_pcffile function seems unnecessary to me, can we instead just take the information from the xyzq? 
        pcf = self.read_pcffile()
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][k]) / 0.52917721)
                curr_mod_charge = (
                    float(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][3])) - m2charges[i][j]
                )
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    distance = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / distance
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    distance = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / distance
        # now linkcorr. linkcorr are last m2*2 entries in pcf
        m2count = 0
        linkstart = len(pcf) - 2 * len(list(_flatten(self.system.list_atoms_m2)))
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
                curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
                m2count += 1
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    distance = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / distance
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    distance = np.linalg.norm(v12)
                    linkenergy += z1 * curr_mod_charge / distance
        # now, add the correction of energy for the link atoms. currently only C-C bond cuts supported.
        for i in range(0, len(self.system.list_coordinates_linkatoms)):
            v1 = [
                self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
            ]
            _flattened = list(_flatten(self.system.list_atoms_q1))
            v2 = [
                self.system.array_xyzq_current[int(_flattened[i]) - 1][0] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][1] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][2] / 0.52917721,
            ]
            v12 = np.array(v2) - np.array(v1)
            distance = np.linalg.norm(v12)
        distance = distance * 0.7409471631
        energycorr = self.databasecorrection("ENERGY", "aminoacid_CACB", distance)

        # sign inverted due to correction convention (subtracting)
        linkenergy -= energycorr

        return linkenergy
    
    def read_pcffile(self):
        pcf = []
        with open(self.PCF.pcf_filename) as ifile:
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

    def get_m2charges(self):
        m2charges = []
        for count in range(len(self.system.list_atoms_m1)):
            m2chargeline = []
            for i in range(0, len(self.system.list_atoms_m2[count])):
                m2chargeline.append(float(self.system.array_xyzq_current[int(self.system.list_atoms_m2[count][i]) - 1][3]))
            m2charges.append(m2chargeline)
        return m2charges
        
    def get_qmenergy(self):

        # logger(logfile, "Extracting QM energy.\n")
        qmenergy = 0.0
        qm_corrdata = []
        if str(self.dict_input_userparameters['program']) == "G16":    
            with open(str(self.str_inputfile_qm + ".log")) as ifile:
                for line in ifile:
                    
                    match = []
                    match2 = []
                    match2 = re.search(
                        r"\sTD[=(\s]", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                    )
                    if not match2:
                        match2 = re.search(
                            r"^TD[=(\s]", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                        )
                    if not match2:
                        match2 = re.search(
                            r"\sTD$", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
                        )
                    if not match2:
                        match2 = re.search(
                            r"^TD$", self.dict_input_userparameters['extra'].upper(), flags=re.MULTILINE
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
                        # logger(logfile, "Obtaining charge self-interaction...\n")
                        pcf_self_pot = self.read_pcf_self()
                        # logger(
                        #     logfile, "Done: {:>20.10f} a.u.\n".format(float(pcf_self_pot))
                        # )
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
        # logger(logfile, "QM energy is " + str(float(qmenergy)) + " a.u..\n")
        return qmenergy, qm_corrdata
    
    def read_pcf_self(self):
        pcf_self = 0.0
        with open(self.str_inputfile_qm + ".log") as ifile:
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
    
    
    def get_mmenergy(self):
        prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        mmenergy = 0.0
        # logger(logfile, "Extracting MM energy.\n")
        # p = subprocess.Popen(
        #     [
        #         prefix,
        #         "energy",
        #         "-f",
        #         edrname,
        #         "-o",
        #         str(edrname + ".xvg"),
        #         "-backup",
        #         "no",
        #     ],
        #     stdout=subprocess.PIPE,
        #     stdin=subprocess.PIPE,
        #     stderr=subprocess.STDOUT,
        # )
        # p.communicate(input=b"11\n\n")
        
        with open(str(self.edrname + ".xvg")) as ifile:
            for line in ifile:
                match = re.search(
                    r"^    0.000000\s*([-]*\d+.\d+)\n", line, flags=re.MULTILINE
                )
                if match:
                    mmenergy = float(match.group(1)) * 0.00038087988
                    break
        # logger(logfile, "MM energy is " + str(float(mmenergy)) + " a.u..\n")
        return mmenergy

    def read_forces(self):
        qmforces = []
        mmforces = []
        qmforces = self.get_qmforces_au()
        # logger(logfile, str("QM forces read.\n"))
        mmforces = self.get_mmforces_au()
        # logger(logfile, str("MM forces read.\n"))
        if self.system.linkcorrlist:
            linkcorrforces = self.get_linkforces_au()
        else:
            linkcorrforces = 0.0
        # logger(logfile, str("Forces for link atom correction read.\n"))
        total_force = np.array(qmforces) + np.array(mmforces) - np.array(linkcorrforces)
        # logger(logfile, str("Total forces obtained.\n"))
        if (self.dict_input_userparameters['jobtype'] != "SINGLEPOINT") and (len(self.system.list_atoms_active) != 0):
            total_force=self.remove_inactive(total_force) #in SP case I don't have inactive atoms
            # logger(logfile, str("Deleted forces of inactive atoms.\n"))
        return total_force
    
    def get_qmforces_au(self):

        qmforces = []
        qmonlyforcelist = []
        pcf_grad = []

        if self.dict_input_userparameters['program'] == "G16":
            insert = ""
            if (int(self.system.int_step_current) != 0): 
                insert = str("." + str(self.system.int_step_current))  
            # logger(logfile,"Reading QM forces using file: "+str(jobname + insert + ".gjf.log")+" and "+ str(jobname + insert + ".fort.7")+"\n")
            qmlogfile = str(self.dict_input_userparameters['jobname'] + insert + ".gjf.log")
            fortfile = str(self.dict_input_userparameters['jobname'] + insert + ".fort.7")

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
        with open(self.class_topology_qmmm.qmmm_topology) as ifile:
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
                                        not in np.array(self.system.list_atoms_qm).astype(int)
                                    )
                                    and (int(match.group(1)) not in np.array(self.system.list_atoms_m1).astype(int))
                                ):
                                    curr_charge = float(match.group(7))
                                    qmforces.append(
                                        [
                                            pcf_grad[count][0] * curr_charge,
                                            pcf_grad[count][1] * curr_charge,
                                            pcf_grad[count][2] * curr_charge,
                                        ]
                                    )
                                    count += 1
                                elif match and int(match.group(1)) in np.array(
                                    self.system.list_atoms_qm
                                ).astype(int):
                                    qmforces.append(qmonlyforcelist[qmcount])
                                    qmcount += 1
                                elif match and int(match.group(1)) in np.array(self.system.list_atoms_m1).astype(
                                    int
                                ):
                                    qmforces.append(
                                        qmonlyforcelist[m1count + len(self.system.list_atoms_qm)]
                                    )
                                    m1count += 1
                                match = re.search(r"^\s*\n", line, flags=re.MULTILINE)
                                if match:
                                    break
                            break
                    break
        return qmforces


    def get_mmforces_au(self):
        prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        mmforces = []
        insert = ""
        if int(self.system.int_step_current) != 0:
            insert = str("." + str(self.system.int_step_current))
        # logger(logfile,"Reading MM forces using file: "+str(self.dict_input_userparameters['jobname'] + insert + ".trr/.tpr/.tpr")+"\n")
        trrname = str(self.dict_input_userparameters['jobname'] + insert + ".trr")
        tprname = str(self.dict_input_userparameters['jobname'] + insert + ".tpr")
        xvgname = str(self.dict_input_userparameters['jobname'] + insert + ".xvg")
        # p = subprocess.Popen(
        #     [
        #         prefix,
        #         "traj",
        #         "-fp",
        #         "-f",
        #         trrname,
        #         "-s",
        #         tprname,
        #         "-of",
        #         xvgname,
        #         "-xvg",
        #         "none",
        #         "-backup",
        #         "no",
        #     ],
        #     stdout=subprocess.PIPE,
        #     stdin=subprocess.PIPE,
        #     stderr=subprocess.STDOUT,
        # )
        # p.communicate(input=b"0\n")

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


    def get_linkforces_au(self):

        linkforces = []
        # Force Coulomb: z1*z2*(distance along coord)/(distance between charges)**3
        for element in self.system.array_xyzq_current:  # this is just to count an entry for each atom!
            linkforces.append([0.0, 0.0, 0.0])
        m2charges = self.get_m2charges()
        for element in self.system.linkcorrlist:
            z1 = 0.0
            v1 = []
            v2 = []
            if int(element[0]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[0]):
                            z1 = float(m2charges[i][j])
                            v1 = [
                                self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z1 != 0.0:
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[0]):
                        z1 = float(self.qm_corrdata[i][2])
                        v1 = [
                            self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[0]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[0]):
                        z1 = float(self.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v1 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z1 = float(self.system.array_xyzq_current[int(element[0]) - 1][3])
                v1 = [
                    self.system.array_xyzq_current[int(element[0]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[0]) - 1][2] / 0.52917721,
                ]
            z2 = 0.0
            if int(element[1]) in _flatten(self.system.list_atoms_m2):
                for i in range(0, len(self.system.list_atoms_m2)):
                    for j in range(0, len(self.system.list_atoms_m2[i])):
                        if int(self.system.list_atoms_m2[i][j]) == int(element[1]):
                            z2 = float(m2charges[i][j])
                            v2 = [
                                self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                                self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                            ]
                            break
                    if z2 != 0.0:
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_qm).astype(int):
                for i in range(0, len(self.system.list_atoms_qm)):
                    if int(self.system.list_atoms_qm[i]) == int(element[1]):
                        z2 = float(self.qm_corrdata[i][2])
                        v2 = [
                            self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                            self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
                        ]
                        break
            elif int(element[1]) in np.array(self.system.list_atoms_m1).astype(int):
                for i in range(0, len(self.system.list_atoms_m1)):
                    if int(self.system.list_atoms_m1[i]) == int(element[1]):
                        z2 = float(self.qm_corrdata[i + len(self.system.list_atoms_qm)][2])
                        v2 = [
                            self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                            self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
                        ]
                        break
            else:
                z2 = float(self.system.array_xyzq_current[int(element[1]) - 1][3])
                v2 = [
                    self.system.array_xyzq_current[int(element[1]) - 1][0] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][1] / 0.52917721,
                    self.system.array_xyzq_current[int(element[1]) - 1][2] / 0.52917721,
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
        pcf = self.read_pcffile()
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i])):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][k]) / 0.52917721)
                curr_mod_charge = (
                    float(float(pcf[int(self.system.list_atoms_m2[i][j]) - 1][3])) - m2charges[i][j]
                )
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_qm[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                        linkforces[int(self.system.list_atoms_m2[i][j]) - 1][l] -= (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_m1[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                        linkforces[int(self.system.list_atoms_m2[i][j]) - 1][l] -= (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
        m2count = 0
        linkstart = len(pcf) - 2 * len(list(_flatten(self.system.list_atoms_m2)))
        for i in range(0, len(self.system.list_atoms_m2)):
            for j in range(0, len(self.system.list_atoms_m2[i]) * 2):
                curr_mod = []
                for k in range(0, 3):
                    curr_mod.append(float(pcf[int(linkstart) + m2count][k]) / 0.52917721)
                curr_mod_charge = float(float(pcf[int(linkstart) + m2count][3]))
                m2count += 1
                for k in range(0, len(self.system.list_atoms_qm)):
                    v1 = [
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][0] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][1] / 0.52917721,
                        self.system.array_xyzq_current[int(self.system.list_atoms_qm[k]) - 1][2] / 0.52917721,
                    ]
                    z1 = float(self.qm_corrdata[k][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_qm[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
                for k in range(0, len(self.system.list_coordinates_linkatoms)):
                    v1 = [
                        self.system.list_coordinates_linkatoms[k][0] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][1] / 0.52917721,
                        self.system.list_coordinates_linkatoms[k][2] / 0.52917721,
                    ]
                    z1 = float(self.qm_corrdata[k + len(self.system.list_atoms_qm)][2])
                    v12 = np.array(v1) - np.array(curr_mod)
                    dist = np.linalg.norm(v12)
                    for l in range(0, 3):
                        linkforces[int(self.system.list_atoms_m1[k]) - 1][l] += (
                            z1 * curr_mod_charge * v12[l] / (dist * dist * dist)
                        )
        for i in range(0, len(self.system.list_coordinates_linkatoms)):
            v1 = [
                self.system.list_coordinates_linkatoms[i][0] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][1] / 0.52917721,
                self.system.list_coordinates_linkatoms[i][2] / 0.52917721,
            ]
            _flattened = list(_flatten(self.system.list_atoms_q1))
            v2 = [
                self.system.array_xyzq_current[int(_flattened[i]) - 1][0] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][1] / 0.52917721,
                self.system.array_xyzq_current[int(_flattened[i]) - 1][2] / 0.52917721,
            ]
            v12 = np.array(v2) - np.array(v1)
            dist = np.linalg.norm(v12) / 0.71290813568205
            unit_vector_v12 = v12 / np.linalg.norm(v12)
            dist = dist * 0.5282272551
            forcecorr = self.databasecorrection("FORCES", "aminoacid_CACB", dist)
            for j in range(0, 3):
                linkforces[int(_flattened[i]) - 1][j] += -unit_vector_v12[j] * forcecorr * 0.5
                linkforces[int(self.system.list_atoms_m1[i]) - 1][j] += unit_vector_v12[j] * forcecorr * 0.5
        return linkforces


    def make_qm_input(self):
        if self.dict_input_userparameters['g16cmd'] == 'g16':
            insert = ""
            oldinsert = ""
            if int(self.system.int_step_current) > 0:
                insert = str("." + str(int(self.system.int_step_current) ))
                if int(self.system.int_step_current) > 1:
                    oldinsert = str("." + str(int(self.system.int_step_current) - 1))
            gaufile = str(self.dict_input_userparameters['jobname'] + insert + ".gjf")
            chkfile = str(self.dict_input_userparameters['jobname'] + insert + ".chk")
            oldchkfile = str(self.dict_input_userparameters['jobname'] + oldinsert + ".chk")

            # XX AJ I have no idea about NMA stuff, is it sufficient here to check for jobtype NMA or can we have the NMA flag also for other jobtypes??
            if self.dict_input_userparameters['jobtype'] == 1:
                oldchkfile = str(self.dict_input_userparameters['jobname'] + ".chk")

            #   Filter xyzq For Coordinates (No Charges)
            fullcoords = filter_xyzq(self.system.array_xyzq_current, list(range(1,self.system.int_number_atoms)), charges=False)

            # XX AJ I don't understand this stdout stuff, Florian could you check if that's correct?
            original_stdout = sys.stdout

            # XX AJ what is the benefit of using stdout instead of f.write?
            with open('top_info.txt', 'w') as f:
                sys.stdout = f # Change the standard output to the file we created.
                print('QMMMTOPINFO')
                for i in self.system.list_atom_elements:
                    print(str(i))
                sys.stdout = original_stdout # Reset the standard output to its original value
            with open('gro_info.txt', 'w') as f:
                sys.stdout = f # Change the standard output to the file we created.
                print('GROINFO')
                for i in fullcoords:
                    print(str(i))
                sys.stdout = original_stdout # Reset the standard output to its original value

            # XX AJ change this file writing if we use default parameters or so
            with open(gaufile, "w") as ofile:
                ofile.write("%NPROCSHARED=" + str(self.dict_input_userparameters['cores']) + "\n")
                ofile.write("%MEM=" + str(self.dict_input_userparameters['memory']) + "MB\n")
                ofile.write("%CHK=" + chkfile + "\n")
                if int(self.system.int_step_current) != 0 or self.nmaflag == 1:
                    ofile.write("%OLDCHK=" + oldchkfile + "\n")
                ofile.write("#P " + str(self.dict_input_userparameters['method']))
                if str(self.dict_input_userparameters['basis']) != "NONE":
                    ofile.write("/" + str(self.dict_input_userparameters['basis']))
                if str(self.dict_input_userparameters['extra']) != "NONE":
                    ofile.write(" " + str(self.dict_input_userparameters['extra']))
                if int(self.system.int_step_current) != 0 or self.nmaflag == 1:
                    ofile.write(" guess=read")
                ofile.write(
                    " nosymm gfinput gfprint force charge guess=huckel punch=derivatives iop(3/33=1,6/7=3) prop(field,read) pop=esp\n"
                )
                ofile.write(
                    "\nQM part of step "+str(self.system.int_step_current)+"\n\n"
                    + str(int(self.dict_input_userparameters['charge']))
                    + " "
                    + str(int(self.dict_input_userparameters['multiplicity']))
                    + "\n"
                )
                count = 0
                for element in fullcoords:
                    if int(count + 1) in np.array(self.system.list_atoms_qm).astype(int):
                        ofile.write(
                            "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                                str(self.system.list_atom_elements[count]),
                                float(element[0]),
                                float(element[1]),
                                float(element[2]),
                            )
                        )
                    count += 1
                for element in self.system.list_coordinates_linkatoms:
                    ofile.write(
                        "{:<2s} {:>12.6f} {:>12.6f} {:>12.6f}\n".format(
                            str("H"), float(element[0]), float(element[1]), float(element[2])
                        )
                    )
                ofile.write("\n")


                # XX AJ we're looping through the pcffile twice here, once for coordinates and charges and once for coordinates where we want the derivate, reduce to one loop later
                with open(self.PCF.pcf_filename) as pcffile:
                    for line in pcffile:
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

                with open(self.PCF.pcf_filename) as pcffile:
                    for line in pcffile:
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

    def run_g16(self):
        insert = ""
        if int(self.system.int_step_current) > 0:
            insert = str("." + str(int(self.system.int_step_current)))

        if not os.path.isfile(str(self.str_inputfile_qm) + ".log"):
            # logger(logfile, "Running G16 file.\n")
            # XX AJ commented out until testing 
            # subprocess.call([g16cmd, str(qmfile)])
            logname = self.str_inputfile_qm[:-3]
            logname += "log"
            # os.rename(logname, str(self.dict_input_userparameters['jobname'] + insert + ".gjf.log"))
            os.rename("fort.7", str(self.dict_input_userparameters['jobname'] + insert + ".fort.7"))
            # logger(logfile, "G16 Done.\n")
        else:
            # logger(
            #     logfile,
            #     "NOTE: Using existing G16 files, skipping calculation for this step.\n",
            # )
            print('XX')
        if not os.path.isfile(self.dict_input_userparameters['jobname'] + insert + ".fort.7"):
            if not os.path.isfile("fort.7"):
                # logger(
                #     logfile,
                #     "No fort.7 file was created by the last Gaussian run! Exiting.\n",
                # )
                # exit(1)
                print('XX error')
            os.rename("fort.7", str(self.dict_input_userparameters['jobname'] + insert + ".fort.7"))
            # logger(
            #     logfile,
            #     "WARNING: Had to rename fort.7 file but not the log file. MAKE SURE THAT THE FORT.7 FILE FITS TO THE LOG FILE!\n",
            # )

    def make_gmx_inp(self):
        self.prefix =  self.dict_input_userparameters['gmxpath'] + self.dict_input_userparameters['gmxcmd']

        insert = ""
        if int(self.system.int_step_current) > 0:
            insert = str("." + str(int(self.system.int_step_current)))

        self.mdpname = str(self.dict_input_userparameters['jobname'] + ".mdp")
        self.groname = str(self.dict_input_userparameters['jobname'] + ".boxlarge.g96")
        self.ndxname = str(self.class_topology_qmmm.qmmm_topology + ".ndx")
        self.tprname = str(self.dict_input_userparameters['jobname'] + insert + ".tpr")

        #   Calculate The Maximum Eucledian Distance Between Any Two Atoms
        array_coordinates_all = self.system.array_xyzq_current[:,:3]
        self.float_distance_max = np.max(np.linalg.norm(array_coordinates_all[np.newaxis, :, :] - array_coordinates_all[:, np.newaxis, :], axis=-1))


        self.write_mdp()
        self.update_gro_box()

        # XX AJ commented out until testing 
        # subprocess.call(
        #     [
        #         prefix,
        #         "grompp",
        #         "-p",
        #         str(self.class_topology_qmmm.qmmm_topology),
        #         "-c",
        #         str(self.groname),
        #         "-n",
        #         str(self.ndxname),
        #         "-f",
        #         str(self.mdpname),
        #         "-o",
        #         str(self.tprname),
        #         "-backup",
        #         "no",
        #     ]
        # )
        # subprocess.call(["rm", "mdout.mdp"])


    def write_mdp(self):
            if self.dict_input_userparameters['rcoulomb'] == 0:
                self.dict_input_userparameters['rcoulomb'] = self.float_distance_max
            if self.dict_input_userparameters['rvdw'] == 0:
                self.dict_input_userparameters['rvdw'] = self.dict_input_userparameters['rcoulomb']

            with open(self.mdpname, "w") as ofile:
                ofile.write(
                    "title               =  Yo\ncpp                 =  /usr/bin/cpp\nconstraints         =  none\nintegrator          =  md\ndt                  =  0.001 ; ps !\nnsteps              =  1\nnstcomm             =  0\nnstxout             =  1\nnstvout             =  1\nnstfout             =   1\nnstlog              =  1\nnstenergy           =  1\nnstlist             =  1\nns_type             =  grid\nrlist               =  "
                )
                ofile.write(str(float(self.dict_input_userparameters['rcoulomb'])))
                ofile.write(
                    "\ncutoff-scheme = group\ncoulombtype    =  cut-off\nrcoulomb            =  "
                )
                ofile.write(str(float(self.dict_input_userparameters['rcoulomb'])))
                ofile.write("\nrvdw                =  ")
                ofile.write(str(float(self.dict_input_userparameters['rvdw'])))
                if self.dict_input_userparameters['useinnerouter']:
                    ofile.write(
                    "\nTcoupl              =  no\nfreezegrps          =  OUTER\nfreezedim           =  Y Y Y\nenergygrps          =  QM INNER OUTER\nenergygrp-excl = QM QM INNER OUTER OUTER OUTER\nPcoupl              =  no\ngen_vel             =  no\n"
                    )
                else:
                    ofile.write(
                    "\nTcoupl              =  no\nenergygrps          =  QM\nenergygrp-excl = QM QM\nPcoupl              =  no\ngen_vel             =  no\n"
                    )

    def update_gro_box(self):
        with open(self.groname, "w") as ofile:
            with open(self.dict_input_userparameters['coordinatefile']) as ifile:
                # logger(
                #     logfile,
                #     str(
                #         "Finding a larger .gro box size to avoid problems with .mdp input...\n"
                #     ),
                # )
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
                                # logger(
                                #     logfile,
                                #     "\n\nError: In "
                                #     + str(gro)
                                #     + " box vectors were expected but not found. Exiting. Line was:\n",
                                # )
                                # logger(logfile, line)
                                exit(1)
                            else:
                                bv = [
                                    float(match.group(1)) + 10.0 * self.float_distance_max,
                                    float(match.group(2)) + 10.0 * self.float_distance_max,
                                    float(match.group(3)) + 10.0 * self.float_distance_max,
                                ]
                                ofile.write(
                                    " {:>15.9f} {:>15.9f} {:>15.9f}\nEND\n".format(
                                        float(bv[0]), float(bv[1]), float(bv[2])
                                    )
                                )
                            break
                        break
        # logger(logfile, str("Done.\n"))

    def run_gmx(self):

        insert = ""
        if int(self.system.int_step_current) != 0:
            insert = str("." + str(int(self.system.int_step_current)))

        # logger(logfile, "Running Gromacs file.\n")
        self.trrname = str(self.dict_input_userparameters['jobname'] + insert + ".trr")
        self.xtcname = str(self.dict_input_userparameters['jobname'] + insert + ".xtc")
        self.outname = str(self.dict_input_userparameters['jobname'] + insert + ".out.gro")
        self.gmxlogname = str(self.dict_input_userparameters['jobname'] + insert + ".gmx.log")
        self.edrname = str(self.dict_input_userparameters['jobname'] + insert + ".edr")

        # XX AJ commented out until testing 
        # subprocess.call(
        # [
        #     self.prefix,
        #     "mdrun",
        #     "-s",
        #     self.tprname,
        #     "-o",
        #     self.trrname,
        #     "-c",
        #     self.outname,
        #     "-x",
        #     self.xtcname,
        #     "-g",
        #     self.gmxlogname,
        #     "-e",
        #     self.edrname,
        #     "-backup",
        #     "no",
        # ]
        # )
        # # os.remove(outname)

    def databasecorrection(self, energy_or_force, cut, distance):
        
        # XX AJ check what this connection object is
        conn = sqlite3.connect(self.basedir + "/correction_database/database.sqlite")

        # check if method exist in database
        method_set = conn.cursor()
        method_set.execute(
            "SELECT * FROM "
            + cut
            + ' WHERE forcefield="'
            + self.dict_input_userparameters['forcefield']
            + '" AND method="'
            + self.dict_input_userparameters['method']
            + '" AND basisset ="'
            + self.dict_input_userparameters['basis']
            + '"'
        )
        method_set_value = method_set.fetchall()
        if len(method_set_value) == 0:
            cut = "aminoacid_CACB"
            self.dict_input_userparameters['forcefield'] = "amberGS"
            self.dict_input_userparameters['method'] = "CAM-B3LYP"
            self.dict_input_userparameters['basis'] = "6-31++G**"
            # logger(
            #     logfile,
            #     "Unexisted method in correction database, changing to default correction method...\n",
            # )

        c = conn.cursor()
        c.execute(
            "SELECT * FROM "
            + cut
            + ' WHERE forcefield="'
            + self.dict_input_userparameters['forcefield']
            + '" AND method="'
            + self.dict_input_userparameters['method']
            + '" AND basisset ="'
            + self.dict_input_userparameters['basis']
            + '"'
        )
        db_values = c.fetchall()[0]

        conn.close()
        returnvalue = 0
        if len(db_values) > 0:
            if energy_or_force == "ENERGY":
                if self.dict_input_userparameters['databasefit'] == "POLY":
                    returnvalue = (
                        db_values[5] * distance * distance * distance
                        + db_values[6] * distance * distance
                        + db_values[7] * distance
                        + db_values[8]
                    )
                elif self.dict_input_userparameters['databasefit'] == "MORSE":
                    returnvalue = (
                        db_values[9]
                        * (
                            np.exp(-2 * db_values[10] * (distance - db_values[11]))
                            - 2 * np.exp(-db_values[10] * (distance - db_values[11]))
                        )
                        + db_values[12]
                    )
                elif self.dict_input_userparameters['databasefit'] == "NO":
                    returnvalue = 0
                    # logger(logfile, "No energy correction.\n")
            elif energy_or_force == "FORCES":
                if self.dict_input_userparameters['databasefit'] == "POLY":
                    returnvalue = (
                        db_values[13] * distance * distance * distance
                        + db_values[14] * distance * distance
                        + db_values[15] * distance
                        + db_values[16]
                    )
                elif self.dict_input_userparameters['databasefit'] == "MORSE":
                    returnvalue = (
                        db_values[17]
                        * (
                            np.exp(-2 * db_values[18] * (distance - db_values[19]))
                            - 2 * np.exp(-db_values[18] * (distance - db_values[19]))
                        )
                        + db_values[20]
                    )
                elif self.dict_input_userparameters['databasefit'] == "NO":
                    returnvalue = 0
                    # logger(logfile, "No force correction.\n")
        return returnvalue
    
    def write_output(self, energy_file="oenergy.txt", forces_file="oforces.txt"):
        qmenergy, mmenergy, linkcorrenergy, total_energy = self.energies
        
        if self.system.int_step_current == 0:
            file_flag = "w"
        else:
            file_flag = "a+"

        oenergy = open(energy_file, file_flag)
        if self.system.int_step_current == 0:
            oenergy.write("Step\tQM\t\tMM\t\tLink\t\tTotal\n")
            oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (self.system.int_step_current, qmenergy, mmenergy, linkcorrenergy, total_energy))
        else:
            
            oenergy.write("%d\t%f\t%f\t%f\t%f\n" % (self.system.int_step_current, qmenergy, mmenergy, linkcorrenergy, total_energy))
        oenergy.close()

        oforce = open(forces_file, file_flag)
        oforce.write("Step%d\n"%self.system.int_step_current)
        for i in range(len(self.total_force)):
            oforce.write(
                "%d %.8f %.8f %.8f\n"
                % ((i + 1), self.total_force[i][0], self.total_force[i][1], self.total_force[i][2])
            )
        oforce.write("\n")
        oforce.close()