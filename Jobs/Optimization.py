#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Run Optimizations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-08-20'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os
import math
import sqlite3
import sys
import copy
import subprocess
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries
import Jobs.Singlepoint as SP

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators.GeneratorGeometries import propagate_dispvec
from Generators._helper import filter_xyzq, _flatten

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Optimization():

    '''
    This Class Performs An Optimization
    '''

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base) -> None:
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        XX\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        input_dict: dict -> Dictionary Of The User Parameter Input \\
        class_system: class -> Class Object Of The System \\
        class_topology: class -> Class Object Of The Topology \\
        class_pcf: class -> Class Object Of The Pointchargefield \\
        str_directory_base: str -> Directory Path \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.str_directory_base = str_directory_base

        #   XX AJ check how to deal with nma flag later
        self.nmaflag = 0

        #define done XX AJ evaluate this later
        STEPLIMIT = 0
        FTHRESH = 1
        STEPSIZE = 2

        #   Perform Initial Singlepoint Calculation
        self.singlepoint = SP.Singlepoint(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base)

        '''
        temporary documentation on how variables are called in old and new version
        ---------------------------------------------------------------------
        old                         new
        ---------------------------------------------------------------------
        curr_energy                 self.singlepoint.total_energy
        total_force                 self.singlepoint.total_force
        last_forces                 self.list_forces_all_steps
        '''

        #   Setting Up Variables
        self.list_forces_max_all_steps = []
        self.list_forces_all_steps = []
        self.list_forces_all_steps.append(self.singlepoint.total_force)

        self.list_energies_all_steps = []
        self.list_energies_all_steps.append([self.singlepoint.qmenergy, self.singlepoint.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy])

        # XX AJ I'm not sure if the lists needs 2 elements, probably for scan or BFGS, I will come back to that
        self.list_xyzq_all_steps = [self.system.array_xyzq_current, self.system.array_xyzq_current]

        # XX AJ at this point I don't think we still need this function, but I'll keep these comments for now in case I'm wrong
        # self.force_clean = self.make_clean_force()

        #   Check If Maximum Force Is Below Threshold
        float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
        if abs(float_force_max) < self.dict_input_userparameters['f_thresh']:
            self.bool_opt_done = FTHRESH
        else:
            self.bool_opt_done = False
        self.list_forces_max_all_steps.append(float_force_max)

        #   Starting Optimization Cycles
        while not self.bool_opt_done and self.system.int_step_current <= self.dict_input_userparameters['maxcycle']:
            self.system.int_step_current += 1

            #   Prepare New Input
            self.singlepoint.groname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".g96")
            self.singlepoint.tprname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".tpr")
            self.singlepoint.trrname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".trr")
            self.singlepoint.xtcname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".xtc")
            self.singlepoint.outname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".out.gro")
            self.singlepoint.gmxlogname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gmx.log")
            self.singlepoint.edrname = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr")

            # gaufile = str(jobname + insert + ".gjf")
            # chkfile = str(jobname + insert + ".chk")
            # oldchkfile = str(jobname + oldinsert + ".chk")
            self.PCF.pcf_filename = str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".pointcharges")


            if self.dict_input_userparameters['jobtype'] == "SCAN" :
                pass # XX AJ add scan later
                # dispvec = propagate_dispvec(propagator, xyzq, new_xyzq, total_force, last_forces, stepsize, self.system.int_step_current, True, qmmmInputs.scan_atoms)
            else :
                dispvec = propagate_dispvec(self.dict_input_userparameters['propagator'], self.list_xyzq_all_steps, self.list_forces_all_steps, self.list_forces_max_all_steps[-1], self.dict_input_userparameters['stepsize'], self.system.int_step_current)
            #    write_dispvec(dispvec, curr_step, count_trash) Simon implemented this to check if the dispvec is correct!

            #   Apply Displacement
            list_xyzq_new = self.list_xyzq_all_steps[-1] + np.append(dispvec, np.zeros((len(dispvec),1)), axis=1)
            self.list_xyzq_all_steps.append(list_xyzq_new)

            #   Update Pointchargefiel
            self.PCF.make_pcf()

            #Run SP
            self.singlepoint.run_calculation()

            #   Update Forces And Energies
            self.list_forces_all_steps.append(self.singlepoint.total_force)
            self.list_energies_all_steps.append([self.singlepoint.qmenergy, self.singlepoint.mmenergy, self.singlepoint.linkcorrenergy, self.singlepoint.total_energy])

            #   Check If The Total Energy Improved
            if self.list_energies_all_steps[-1][-1] > self.list_energies_all_steps[-2][-1]:
                self.bool_energy_improved = False
                self.dict_input_userparameters['stepsize'] *= 0.2

                #   Remove Files
                os.remove(self.singlepoint.trrname)
                os.remove(self.singlepoint.tprname)
                os.remove(self.singlepoint.gmxlogname)
                os.remove(self.singlepoint.edrname)
                os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".edr.xvg"))
                os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf.log"))
                os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".fort.7"))
                os.remove(str(self.dict_input_userparameters['jobname'] + "." + str(self.system.int_step_current) + ".gjf"))
                os.remove(self.PCF.pcf_filename)

            else:
                self.bool_energy_improved = True
                self.dict_input_userparameters['stepsize'] *= 1.2

            float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
            self.list_forces_max_all_steps.append(float_force_max)

            if abs(float_force_max) < float(self.dict_input_userparameters['f_thresh']):
                # logger(logfile,"Max force (%f) below threshold (%f) Finishing.\n"%(maxforce,f_thresh))
                self.bool_opt_done = FTHRESH
                break

            if float(self.dict_input_userparameters['stepsize']) < 1e-6: #0.000001 a.u.
                self.bool_opt_done = STEPSIZE
                # logger(
                #     logfile,
                #         ("Step became lower than 0.000001 a.u., optimization is considered done for now. " +
                #          "This is the best we can do unless reaching unacceptable numerical noise levels.\n"),
                # )
                break

            if self.bool_energy_improved:
                #   Remove The xyzq List From Two Steps Ago
                self.list_xyzq_all_steps.pop(0)
                #Update for BFGS
                # xyzq = old_qmmmInputs.xyzq
                # new_xyzq = qmmmInputs.xyzq
                #Store previous
                if self.dict_input_userparameters['jobname'] == "SCAN" :
                    # write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step, energy_file="oenergy_%s.txt"%jobname, forces_file="oforces_%s.txt"%jobname)
                    pass
                else:
                    # write_output(qmmmInputs.energies, qmmmInputs.forces, qmmmInputs.qmmmparams.curr_step)
                    pass
                # gro = qmmmInputs.gro            #SIMON
                # logger(logfile, "Due to the decrease of the energy, the structure "+str(gro)+" will be used from now on.\n")
            else:
                #rejected and use previous
                self.system.int_step_current -= 1

                #   Remove Current Parameters
                self.list_xyzq_all_steps.pop()
                self.list_forces_all_steps.pop()
                self.list_forces_max_all_steps.pop()
                # logger(logfile,"Due to the rejection use maximum force: "+str(maxforce)+"\n")
                # logger(logfile,"Due to the increase of the energy, the used structure remains "+str(gro)+".\n")

            # XX logging optimization done due to ...
            
            pass


