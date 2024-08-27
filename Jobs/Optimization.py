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

    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, directory_base) -> None:
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

        self.dict_input_userparameters = dict_input_userparameters
        self.system = class_system
        self.class_topology_qmmm = class_topology
        self.PCF = class_pcf
        self.directory_base = directory_base

        #   XX AJ check how to deal with nma flag later
        self.nmaflag = 0

        #define done XX AJ evaluate this later 
        STEPLIMIT = 0
        FTHRESH = 1
        STEPSIZE = 2
         
        #   Perform Initial Singlepoint Calculation
        self.singlepoint = SP.Singlepoint(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.directory_base)

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
        self.list_forces_all_steps = []
        self.list_forces_all_steps.append(self.singlepoint.total_force)

        # XX AJ at this point I don't think we still need this function, but I'll keep these comments for now in case I'm wrong
        # self.force_clean = self.make_clean_force()

        #   Check If Maximum Force Is Below Threshold
        float_force_max = max(np.max(self.singlepoint.total_force), np.min(self.singlepoint.total_force), key=abs)
        if abs(float_force_max) < self.dict_input_userparameters['f_thresh']:
            self.bool_opt_done = FTHRESH
        else:
            self.bool_opt_done = False

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
                dispvec = propagate_dispvec(self.dict_input_userparameters['propagator'], self.system.array_xyzq_current, self.system.array_xyzq_current, self.list_forces_all_steps, float_force_max, self.dict_input_userparameters['stepsize'], self.system.int_step_current)
            #    write_dispvec(dispvec, curr_step, count_trash) Simon implemented this to check if the dispvec is correct!  
            
            #   Apply Displacement
            self.system.array_xyzq_previous = np.copy(self.system.array_xyzq_current)
            self.system.array_xyzq_current += np.append(dispvec, np.zeros((len(dispvec),1)), axis=1)
                      
            #   Update Pointchargefiel
            self.PCF.make_pcf()

            #Run SP
            self.singlepoint.run_calculation()
            pass


