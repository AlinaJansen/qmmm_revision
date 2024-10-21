#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Run Scan Calculations"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-10-02'

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
from Logging.WriteOutput import Output
from Jobs import QM_job, MM_job
from Jobs import Optimization as OPT
from Generators._helper import filter_xyzq, _flatten
from Generators.GeneratorGeometries import read_gmx_structure_header, read_gmx_structure_atoms, read_gmx_box_vectors, write_g96
from Generators.GeneratorEnergies import GeneratorEnergies, GeneratorForces

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   - Write Output part
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class Scan():

    '''
    This Class Performs A Singlepoint Calculation
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

        self.float_stepsize_scan = float(self.dict_input_userparameters['scanstepsize']) # XX AJ remove float later and put into assessment

        #   XX AJ check how to deal with nma flag later
        self.nmaflag = 0

class Scan_bondlength(Scan):
    def __init__(self, dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base) -> None:
        super().__init__(dict_input_userparameters, class_system, class_topology, class_pcf, str_directory_base)

        #   Create A List For The Scan Atoms
        self.list_atoms_scan = [int(x) for x in self.dict_input_userparameters['scanatomslist'].split(',')]

        
        # self.int_scan_step = self.dict_input_userparameters['scanstep'] XX AJ I think I can delete that, check when I'm done with the scan function
        self.str_directory_scan = f'scanR/R{self.list_atoms_scan[0]}-{self.list_atoms_scan[0]}'
        
        #   Create Scan Directories XX AJ we have to decide if we want to raise errors for existing directories (exist_ok=False) or not (exist_ok=True)
        # os.makedirs(self.str_directory_scan)
        os.makedirs(self.str_directory_scan, exist_ok=True)
        
        for self.int_scan_step in range(int(self.dict_input_userparameters['scanstep']), int(self.dict_input_userparameters['scanstepstotal'])):
            self.str_directory_scan_step = f'{self.str_directory_scan}/step{self.int_scan_step}'
            os.makedirs(self.str_directory_scan_step, exist_ok=True)

            #   Calculate The Displacement
            self.generate_initial_displacement()

            #   Apply The Displacement
            self.system.array_xyzq_initial += self.array_displacement_initial

            #   Update The Grofile
            self._str_coordinatefile_new = f'scanR{int(self.list_atoms_scan[0])}-{int(self.list_atoms_scan[1])}_{self.int_scan_step}.g96'
            self.generate_g96_file_scan()
            self.dict_input_userparameters['coordinatefile'] = self._str_coordinatefile_new

            #   Update The Jobname 
            self.dict_input_userparameters['jobname'] = f'scanR{int(self.list_atoms_scan[0])}-{int(self.list_atoms_scan[1])}_{self.int_scan_step}'

            #   Run Optimization
            OPT_job = OPT.Optimization(self.dict_input_userparameters, self.system, self.class_topology_qmmm, self.PCF, self.str_directory_base)

            #   Write Output
            # XX AJ missing (:
            


    def generate_initial_displacement(self):
        self.array_displacement_initial = np.zeros((len(self.system.array_xyzq_initial),3))
        coords = np.array(self.system.array_xyzq_initial)[:, 0:3]
        array_unit_vector_scan_atoms = (np.array(coords[int(self.list_atoms_scan[1]-1)]) - np.array(coords[int(self.list_atoms_scan[0]-1)])) / np.linalg.norm(np.array(coords[int(self.list_atoms_scan[0]-1)]) - np.array(coords[int(self.list_atoms_scan[1]-1)]))
        line_vec = array_unit_vector_scan_atoms * self.float_stepsize_scan  # ba * stepsize 0.052917721 # 1 b,a.u. = 0.529177249 A 
        self.array_displacement_initial[int(self.list_atoms_scan[1]-1)] += line_vec
        self.array_displacement_initial[int(self.list_atoms_scan[1]-1)] *= 0.052917721

    def generate_g96_file_scan(self):
        #   XX AJ maybe that would be nicer if GeneratorGeometries.py would also be a class ...
        header = read_gmx_structure_header(self.dict_input_userparameters['coordinatefile'])
        list_atom_information = read_gmx_structure_atoms(self.dict_input_userparameters['coordinatefile'])
        list_vectors_box = read_gmx_box_vectors(self.dict_input_userparameters['coordinatefile'])
        write_g96(self._str_coordinatefile_new, header, list_atom_information, self.system.array_xyzq_initial, list_vectors_box)

        
        
        
        
        pass