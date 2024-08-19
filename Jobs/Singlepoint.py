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
import sys
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries
import Generators.GeneratorGeometries as geometry

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators._helper import filter_xyzq

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Singlepoint():

    '''
    This Class Performs A Singlepoint Calculation
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, basedir) -> None:
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
        self.qmenergy, self.mmenergy, self.linkcorrenergy, self.qm_corrdata = self.get_energy(qmfile, edrname, qmmmInputs)
        self.total_energy = self.qmenergy + self.mmenergy - self.linkcorrenergy
        self.energies = (self.qmenergy, self.mmenergy, self.linkcorrenergy, self.total_energy)

        #   Read Forces
        self.total_force = self.read_forces(qm_corrdata,qmmmInputs)

        

    def get_energy(self, qmfile, edrname, qmmmInputs):
        logfile = qmmmInputs.logfile
        qmenergy, qm_corrdata = self.get_qmenergy(qmfile, qmmmInputs)
        mmenergy = self.get_mmenergy(str(edrname), qmmmInputs)
        if qmmmInputs.linkcorrlist:
            linkcorrenergy = self.get_linkenergy_au(qm_corrdata, qmmmInputs)
        else:
            linkcorrenergy = 0.0
        basis = qmmmInputs.qmparams.basis
        methodstring = str(qmmmInputs.qmparams.method)
        if basis != "NONE":
            methodstring += str("/" + str(basis))

        return qmenergy, mmenergy, linkcorrenergy, qm_corrdata

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
                with open(self.system.pcffile) as pcffile:
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

                with open(self.system.pcffile) as pcffile:
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


    def run_g16(self):
        insert = ""
        if int(self.system.int_step_current) > 0:
            insert = str("." + str(int(self.system.int_step_current)))

        if not os.path.isfile(str(self.str_inputfile_qm) + ".log"):
            # logger(logfile, "Running G16 file.\n")
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

        subprocess.call(
        [
            self.prefix,
            "mdrun",
            "-s",
            self.tprname,
            "-o",
            self.trrname,
            "-c",
            self.outname,
            "-x",
            self.xtcname,
            "-g",
            self.gmxlogname,
            "-e",
            self.edrname,
            "-backup",
            "no",
        ]
        )
        # os.remove(outname)