#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-03-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os 
import math
import json
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries
import Generators.GeneratorGeometries as geometry

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators import generate_pcf_from_top as make_pcf
from Generators._helper import _flatten, get_linkatoms_angstrom

#   // TODOS & NOTES //
#   TODO:
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class SystemInfo():

    '''
    This Class Reads And Stores Information About The System
    '''

    def __init__(self, dict_input_userparameters) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Calls Functions To Read Information About The System\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        input_dict: dict -> Dictionary Of The User Parameter Input \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        self.dict_input_userparameters = dict_input_userparameters

        #   Make a list of all topology files
        self.list_topology = self.get_list_topologies(self.dict_input_userparameters['topologyfile'])

        # Read The Different Types Of Molecules In The System
        self.list_molecules = self.read_molecules()

        # Read The Charges For All Atoms
        self.list_charges = []
        for molecule in self.list_molecules:
            self.list_charges.extend(self.readcharges(molecule))

        #   Read All Atom Lists
        self.list_atoms_qm = self.read_list_atoms_input(self.dict_input_userparameters['qmatomslist'])
        self.list_atoms_active = self.read_list_atoms_input(self.dict_input_userparameters['activeatomslist'])
        # XX AJ the inner/outer part will be unnecessary when Nicola adds his part
        if dict_input_userparameters['useinnerouter']:
            self.list_atoms_inner = self.read_list_atoms_input(self.dict_input_userparameters['inneratomslist'])
            self.list_atoms_outer = self.read_list_atoms_input(self.dict_input_userparameters['outeratomslist'])
        #   If we're not using inner/outer we're keeping the lists empty (currently necessary for some functions but that might change later) XX AJ
        else:
            self.list_atoms_inner = []
            self.list_atoms_outer = []

        #   Read Connectivity
        self.list_connectivity_topology = self.read_connectivity_topology(self.dict_input_userparameters['topologyfile'])

        #   Read Initial Geometry
        #   XX AJ I would prefer only one function here independent of the file type and only make that distinction within the function. I'll get back to that later when I'm writing GeneratorGeometries.py
        if self.dict_input_userparameters['coordinatefile'][-4:] == ".gro":
            #logger(logfile, "Reading geometry (.gro)...\n")
            # XX AJ if we rewrite gro files to g96 files, the following function can be deleted and we can read the g96 file with geometry.readg96 afterwards
            self.list_geometry_initial = geometry.readgeo(self.dict_input_userparameters['coordinatefile'])
            # Writing high-precision coordinate file
            # logger(logfile, "Writing high-precision coordinate file...")
            self.write_file_gromacs_highprec(self.dict_input_userparameters['coordinatefile']) # XX temp removed logfile until logfile decision of Florian AJ
            self.dict_input_userparameters['coordinatefile'] = self.dict_input_userparameters['jobname'] + ".g96"
            # logger(logfile, "Done.\n")
        elif self.dict_input_userparameters['coordinatefile'][-4:] == ".g96":
            #logger(logfile, "Reading geometry (.g96)...\n")
            self.list_geometry_initial = geometry.readg96(self.dict_input_userparameters['coordinatefile'])

        self.int_number_atoms = int(len(self.list_geometry_initial)/3)

        #   Create xyzq (Coordinates And Charges) For The Whole System
        #   XX AJ this will also change with Nicolas part
        if dict_input_userparameters['useinnerouter']:
            self.array_xyzq = geometry.make_xyzq_io(self.list_geometry_initial, self.list_charges, self.list_atoms_outer)
        else:
            self.array_xyzq = geometry.make_xyzq(self.list_geometry_initial, self.list_charges)

        #    Read Linkatoms And Next Order Atoms In MM Region; For Further Information On Atom Naming See Documentation
        self.get_list_atoms_m1()
        self.get_list_atoms_m2()

        #   Read coordinates of linkatoms in angstrom
        self.list_coordinates_linkatoms = get_linkatoms_angstrom(self.array_xyzq, self.list_atoms_qm, self.list_atoms_m1, self.list_connectivity_topology, [])
        self.get_list_atoms_link()


    def read_list_atoms_input(self, file_input_atoms) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------
        Reads Atom Indices
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        file_input_atoms: str -> Index File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_atoms_indices: list -> List Of Atom Indices \\
        ------------------------------ \\
        '''

        list_atoms_indices = []
        with open(file_input_atoms, 'r') as atom_file:
            for line in atom_file:

                # Skip Index Group Name
                if "[" in line or "]" in line:
                    continue
                
                list_atom_indices_line = re.findall("\d+", line)
                if list_atom_indices_line:
                    list_atoms_indices.extend(map(int, list_atom_indices_line))
        list_atoms_indices = sorted(np.array(list_atoms_indices).astype(int))
        return list_atoms_indices

    def read_molecules(self) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads List Of Molecules From The Topology File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_molecule: list -> List Of Molecules And Their Amount In The System \\
        ------------------------------ \\
        '''

        list_molecule = []
        with open(self.dict_input_userparameters['topologyfile'], 'r') as file_input:
            bool_match = False
            for line in file_input:
                match = re.search(r"^\[ molecules \]", line, flags=re.MULTILINE)
                if match:
                    bool_match = True
                    break
            if not bool_match:
                #   XX AJ turn into Logging
                print('No "molecules" entry in ' + str(self.top) + " found. Exiting.")
                exit(1)
            for line in file_input:
                # Skip All Comment Lines
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                # Find All Molecule Types And Their Amount
                else:
                    # Extract Lines With A String (Molecule Type, Group 1) And A Number (Molecule Amount, Group 2) Seperated By Whitespace
                    match = re.search(r"^(\S+)\s+(\d+)", line, flags=re.MULTILINE)
                    if match:
                        list_molecule.append([match.group(1), match.group(2)])
                    else:
                        #   XX AJ turn into Logging
                        print("Found an incomprehensible line in molecules list. Exiting.")
                        print("Last line was:")
                        print(line)
                        exit(1)
        return list_molecule

    def readcharges(self, list_molecule_entry) -> list:
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads Charges For All Atoms In A Molecule Type \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        list_molecule_entry: list -> List With Molecule Name And Amount Of This Molecule \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_charges: list -> List Of Charges For Every Atom In This Molecule Types \\
        ------------------------------ \\
        '''

        list_charges = []
        str_topology_current = self.dict_input_userparameters['topologyfile']
        str_molecule_name = list_molecule_entry[0]
        int_molecule_amount = int(list_molecule_entry[1])
        bool_molecule_occurence = self.check_occurence_topology_molecule(str_molecule_name, self.dict_input_userparameters['topologyfile'])

        if not bool_molecule_occurence:
            for topology in self.list_topology:
                bool_molecule_occurence = self.check_occurence_topology_molecule(str_molecule_name, topology)
                if bool_molecule_occurence:
                    str_topology_current = topology
                    break
        if not bool_molecule_occurence:
            print("No charges found for " + str(str_molecule_name) + ". Exiting.")
            exit(1)
        with open(str_topology_current) as input_file:
            for line in input_file:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in input_file:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        matchstring = r"^\s*" + re.escape(str_molecule_name)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            bool_molecule_occurence = True

                            for line in input_file:
                                match = re.search(
                                    r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                )
                                if match:
                                    break
                            break
                        else:
                            bool_molecule_occurence = False
                            break
                    if bool_molecule_occurence:
                        break
            for line in input_file:
                match = re.search(r"^\[", line, flags=re.MULTILINE)
                if match:
                    break
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(
                    r"^\s*\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+([-]*\d+[\.]*[\d+]*)*",
                    line,
                    flags=re.MULTILINE,
                )
                if match:
                    list_charges.append(float(match.group(1)))

        list_charges *= int_molecule_amount

        return list_charges

    def check_occurence_topology_molecule(self, str_molecule_name, str_file_topology) -> bool:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Checks If The Molecule (molecule_name) Is Listed In The Topology File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_molecule_name: str -> Name Of The Molecule \\
        str_file_topology: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        bool_occurence_molecule: bool \\
        ------------------------------ \\
        '''

        with open(str_file_topology) as input_file:
            bool_occurence_molecule = False
            for line in input_file:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in input_file:
                        if line.strip() == '' or line.startswith(';'):
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(str_molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                bool_occurence_molecule = True
                            break
                if bool_occurence_molecule:
                    break

        return bool_occurence_molecule

    def read_connectivity_topology(self, str_file_topology) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_file_topology: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_connectivity: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        int_count = 0
        list_connectivity = []
        for molecule in self.list_molecules:
            current_topology = self.get_topology_current(molecule[0], str_file_topology)
            for _ in range(0, int(molecule[1])):
                int_molecule_atoms = self.get_mollength_direct(molecule[0], current_topology)
                list_molecule_connectivity = self.get_list_connectivity(int_count, molecule[0], current_topology)
                list_connectivity += list_molecule_connectivity
                int_count += int(int_molecule_atoms)

        return list_connectivity

    def get_topology_current(self, str_molecule_name, str_file_topology) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_molecule_name: str -> Name Of The Molecule \\
        str_file_topology: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        str_topology_current: str -> Name Of Current Topology File \\
        ------------------------------ \\
        '''

        str_topology_current = str_file_topology
        bool_occurence_molecule = make_pcf.checkformol(str_molecule_name, str_file_topology)

        if not bool_occurence_molecule:
            for element in self.list_topology:
                bool_occurence_molecule = make_pcf.checkformol(str_molecule_name, element)
                if bool_occurence_molecule:
                    str_topology_current = element
                    break
        if not bool_occurence_molecule:
            # XX AJ replace with logger
            print("No charges found for " + str(str_molecule_name) + ". Exiting.")
            exit(1)

        return str_topology_current

    def get_mollength_direct(self, str_molecule_name, str_file_topology) -> list:
        
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Amount Of Atoms Of The Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_molecule_name: str -> Name Of The Molecule \\
        str_file_topology: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        int_molecule_length: int -> Number Of Atoms In A Molecule \\
        ------------------------------ \\
        '''

        int_molecule_length = 0
        with open(str_file_topology) as input_file:
            for line in input_file:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in input_file:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(str_molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in input_file:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in input_file:
                                                match = re.search(
                                                    r"^;", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    continue
                                                else:
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+[-]*\d+[\.]*[\d+]*",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        int_molecule_length = int(match.group(1))
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                            break
                                break
                    break
        return int_molecule_length

    def get_list_connectivity(self, int_offset, str_molecule_name, str_file_topology) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of One Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        int_offset: int -> Offset Of Atom Indices \\
        str_molecule_name: str -> Name Of The Molecule \\
        str_file_topology: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_connectivity: list -> List Of The Connectivity Of One Molecule \\
        ------------------------------ \\
        '''
                
        list_connectivity = []
        with open(str_file_topology) as input_file:
            for line in input_file:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[ moleculetype \]", line, flags=re.MULTILINE)
                if match:
                    for line in input_file:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(str_molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in input_file:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[ bonds \]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in input_file:
                                                match = re.search(
                                                    r"^;", line, flags=re.MULTILINE
                                                )
                                                if match:
                                                    continue
                                                else:
                                                    match = re.search(
                                                        r"^\s*(\d+)\s+(\d+)",
                                                        line,
                                                        flags=re.MULTILINE,
                                                    )
                                                    if match:
                                                        a = int(match.group(1)) + int(
                                                            int_offset
                                                        )
                                                        b = int(match.group(2)) + int(
                                                            int_offset
                                                        )
                                                        if a > b:
                                                            c = b
                                                            b = a
                                                            a = c
                                                        bool_found_connectivity = False
                                                        for element in list_connectivity:
                                                            if int(element[0]) == a:
                                                                if int(b) not in np.array(
                                                                    element
                                                                ).astype(int):
                                                                    element.append(int(b))
                                                                    bool_found_connectivity = True
                                                        if not bool_found_connectivity:
                                                            list_connectivity.append(
                                                                [int(a), int(b)]
                                                            )
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                            break
                                break
                    break
        return list_connectivity

    def write_file_gromacs_highprec(self, gro_file) -> str:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Writes A g96 File From The gro File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        gro_file: str -> Name Of The gro File \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        filename: str -> Name Of The New g96 File \\
        ------------------------------ \\
        '''
        # XX AJ please ignore everything in this function, I will move it to the GeneratorsGeometry module 
        #  XX not tested yet, bc I didn't use a gro file AJ
        filename=str(self.dict_input_userparameters['jobname'] + ".g96")
        with open(filename,"w") as ofile:
                with open(gro_file) as ifile:
                        count=0
                        for line in ifile:
                                count+=1
                                if count==2:
                                        break
                        ofile.write("TITLE\nProtein\nEND\nPOSITION\n")
                        counter=0
                        finalline=""
                        for line in ifile:
                                match=re.search(r'^([\s\d]{5})(.{5})(.{5})([\s\d]{5})\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)\s*([-]*\d+\.\d*)', line,flags=re.MULTILINE)
                                if not match:
                                        # logger(logfile,str("Successfully wrote " +  str(int(counter)) + " atoms to internal precision format file.\n"))
                                        finalline=line
                                        break
                                else:
                                        counter+=1
                                        ofile.write(str(match.group(1))+" "+ str(match.group(2))+" "+str(match.group(3))+" {:>6d} ".format(int(counter))+ "{:>15.9f} {:>15.9f} {:>15.9f}\n".format(float(match.group(5)),float(match.group(6)),float(match.group(7))))
                        ofile.write("END\nBOX\n")
                        match=re.search(r'^\s*(\d+\.*\d*)\s*(\d+\.*\d*)\s*(\d+\.*\d*)', finalline,flags=re.MULTILINE)
                        if not match:
                                # logger(logfile,str("Unexpected line instead of box vectors. Exiting. Last line:\n"))
                                # logger(logfile,line)
                                exit(1)
                        else:
                                ofile.write(str(" {:>15.9f} {:>15.9f} {:>15.9f}\n".format(float(match.group(1)),float(match.group(2)),float(match.group(3)))))
                        ofile.write("END")

        return filename

    def get_list_atoms_m1(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M1 Atoms \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''
        
        self.list_atoms_m1 = []
        for atom_qm in self.list_atoms_qm:
            list_bonds = self.get_list_atoms_m1qm(atom_qm)
            for bond in list_bonds:
                if (int(bond) not in np.array(self.list_atoms_qm).astype(int)) and (
                    int(bond) not in np.array(self.list_atoms_m1).astype(int)
                ):
                    self.list_atoms_m1.append(bond)

    def get_list_atoms_m1qm(self, int_target_qm) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M1 Atoms For A Specific QM Atom \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        int_target_qm: int -> Index Of Specific QM Atom \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        list_partners: list -> List Of M1 Atoms For This Specific QM Atom \\
        ------------------------------ \\
        '''

        list_partners = []
        for topology in self.list_connectivity_topology:
            bool_occurence_atom = False

            for i in range(0, len(topology)):
                if int(topology[i]) == int(int_target_qm):
                    bool_occurence_atom = True
                    break
            if bool_occurence_atom:
                if int(topology[0]) == int(int_target_qm):
                    for i in range(1, len(topology)):
                        if int(topology[i]) not in np.array(list_partners).astype(int):
                            list_partners.append(int(topology[i]))
                else:
                    if int(topology[0]) not in np.array(list_partners).astype(int):
                        list_partners.append(int(topology[0]))

        return list_partners

    def get_list_atoms_m2(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M2 Atoms \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        self.list_atoms_m2 = []
        for element in self.list_atoms_m1:
            list_atoms_m2_m1 = []
            list_bonds = self.get_list_atoms_m1qm(element)
            for entry in list_bonds:
                if int(entry) not in np.array(self.list_atoms_qm).astype(int):
                    list_atoms_m2_m1.append(entry)
            self.list_atoms_m2.append(list_atoms_m2_m1)

    def get_list_atoms_link(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All Q1, Q2, Q3, M3 And M4 Atom Indices \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        list_link_correction = []
        self.list_atoms_m3 = []
        self.list_atoms_m4 = []
        self.list_atoms_q1 = []
        self.list_atoms_q2 = []
        self.list_atoms_q3 = []

        # get q1 and m2
        for element in self.list_atoms_m1:
            q1line = []
            for entry in self.list_connectivity_topology:
                if int(element) in np.array(entry).astype(int):
                    if int(element) != int(entry[0]):
                        if int(entry[0]) in np.array(self.list_atoms_qm).astype(int):
                            q1line.append(int(entry[0]))
                    else:
                        for i in range(1, len(entry)):
                            if int(entry[i]) in np.array(self.list_atoms_qm).astype(int):
                                q1line.append(int(entry[i]))
            self.list_atoms_q1.append(q1line)
        # get q2
        self.list_atoms_q1 = list(_flatten(self.list_atoms_q1))
        for element in self.list_atoms_q1:
            q2line = []
            for conn in self.list_connectivity_topology:
                if int(element) in np.array(conn).astype(int):
                    if (
                        int(element) != int(conn[0])
                        and (int(conn[0]) in np.array(self.list_atoms_qm).astype(int))
                        and (int(conn[0]) not in np.array([int(x) for x in _flatten(self.list_atoms_q1)]))
                    ):
                        q2line.append(int(conn[0]))
                    elif int(element) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if (int(conn[i]) in np.array(self.list_atoms_qm).astype(int)) and (
                                int(conn[i]) not in np.array([int(x) for x in _flatten(self.list_atoms_q1)])
                            ):
                                q2line.append(int(conn[i]))
            self.list_atoms_q2.append(q2line)
        # get q3
        for element in self.list_atoms_q2:
            q3lineline = []
            for entry in element:
                q3line = []
                for conn in self.list_connectivity_topology:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) in np.array(self.list_atoms_qm).astype(int))
                            and int(conn[0]) not in np.array(self.list_atoms_q1).astype(int)
                            and int(conn[0]) not in np.array([int(x) for x in _flatten(self.list_atoms_q1)])
                        ):
                            q3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if (
                                    int(conn[i]) in np.array(self.list_atoms_qm).astype(int)
                                    and int(conn[i]) not in np.array(self.list_atoms_q1).astype(int)
                                    and int(conn[i]) not in np.array([int(x) for x in _flatten(self.list_atoms_q1)])
                                ):
                                    q3line.append(int(conn[i]))
                q3lineline.append(q3line)
            self.list_atoms_q3.append(q3lineline)
        # get m3
        for element in self.list_atoms_m2:
            m3lineline = []
            for entry in element:
                m3line = []
                for conn in self.list_connectivity_topology:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) not in np.array(self.list_atoms_qm).astype(int))
                            and int(conn[0]) not in np.array(self.list_atoms_m1).astype(int)
                        ):
                            m3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if int(conn[i]) not in np.array(self.list_atoms_qm).astype(
                                    int
                                ) and int(conn[i]) not in np.array(self.list_atoms_m1).astype(int):
                                    m3line.append(int(conn[i]))
                m3lineline.append(m3line)
            self.list_atoms_m3.append(m3lineline)
        # get m4
        for element in self.list_atoms_m3:
            m4linelineline = []
            for entry in element:
                m4lineline = []
                for stuff in entry:
                    m4line = []
                    for conn in self.list_connectivity_topology:
                        if int(stuff) in np.array(conn).astype(int):
                            if (
                                int(stuff) != int(conn[0])
                                and (int(conn[0]) not in np.array(self.list_atoms_qm).astype(int))
                                and int(conn[0]) not in np.array(self.list_atoms_m1).astype(int)
                            ):
                                found = False
                                for morestuff in self.list_atoms_m2:
                                    if int(conn[0]) in np.array(morestuff).astype(int):
                                        found = True
                                        break
                                if not found:
                                    m4line.append(int(conn[0]))
                            elif int(stuff) == int(conn[0]):
                                for i in range(1, len(conn)):
                                    if int(conn[i]) not in np.array(self.list_atoms_qm).astype(
                                        int
                                    ) and int(conn[i]) not in np.array(self.list_atoms_m1).astype(int):
                                        found = False
                                        for morestuff in self.list_atoms_m2:
                                            if int(conn[i]) in np.array(morestuff).astype(
                                                int
                                            ):
                                                found = True
                                                break
                                        if not found:
                                            m4line.append(int(conn[i]))
                    m4lineline.append(m4line)
                m4linelineline.append(m4lineline)
            self.list_atoms_m4.append(m4linelineline)
        # set up link atom charge corr pairs: q1-m2, q1-m3, q2-m2, l-m2, l-m3, l-m4 - l are represented as their m1 counterparts!
        # for all these combinations we later have to calculate the coulomb potential and subtract that because they all are 1-4 excluded
        count = 0
        for element in self.list_atoms_m1:
            linkpairline = []
            for entry in self.list_atoms_m2[count]:
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(self.list_atoms_m3[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(self.list_atoms_m4[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            list_link_correction.append(linkpairline)
            count += 1
        count = 0
        for element in self.list_atoms_q1:
            linkpairline = []
            for stuff in self.list_atoms_m2[count]:
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            for stuff in list(_flatten(self.list_atoms_m3[count])):
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            list_link_correction.append(linkpairline)
            count += 1
        count = 0
        for element in self.list_atoms_q2:
            for entry in element:
                linkpairline = []
                for stuff in self.list_atoms_m2[count]:
                    if int(entry) < int(stuff):
                        linkpairline.append([entry, stuff])
                    else:
                        linkpairline.append([stuff, entry])
                list_link_correction.append(linkpairline)
            count += 1
        array_link_correction = np.array(list(_flatten(list_link_correction))).reshape(-1, 2)
        list_link_correction_final = []
        for element in array_link_correction:
            found = False
            for i in range(0, len(list_link_correction_final)):
                if (
                    list_link_correction_final[i][0] == element[0]
                    and list_link_correction_final[i][1] == element[1]
                ):
                    found = True
                    break
            if found:
                continue
            list_link_correction_final.append(element)
        self.list_link_correction = sorted(list_link_correction_final, key=lambda l: l[0])

    
    def get_list_topologies(self, str_file_topology_input) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------  \\
        Recursively Looks For Other Topology Files Referred To In The Original File \\
        And Adds Them To The List Of Topology Files \\
        ------------------------------ \\
        INPUT: \\
        ---------------  \\
        str_file_topology_input: str -> Name Of Topology \\
        ------------------------------ \\
        RETURN: \\
        ---------------  \\
        list_files_topologies: list -> List Of Topology Files \\
        ------------------------------ \\
        '''
        
        list_files_topologies = []
        with open(str_file_topology_input) as input_file:
            for line in input_file:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^#include\s+\"(\S+)\"", line, flags=re.MULTILINE)
                if match:
                    match2 = re.search("ffbonded", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("ffnonbonded", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("forcefield.itp", match.group(1))
                    if match2:
                        continue
                    match2 = re.search("posre.itp", match.group(1))
                    if match2:
                        continue
                    foundname = match.group(1)
                    check = os.path.isfile(foundname)
                    if not check:
                        foundname = os.path.join(self.dict_input_userparameters['gmxtop_path'], *foundname.strip('/').split('/'))
                        check = os.path.isfile(foundname)
                        if not check:
                            print("File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting.")
                            exit(1)
                    list_files_topologies.append(foundname)

                    list_files_topologies.extend(self.get_list_topologies(foundname))
        return list_files_topologies
    
    def get_atoms(self) -> list:

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------  \\
        Reads Elements For All Atoms \\
        ------------------------------ \\
        INPUT: \\
        ---------------  \\
        None \\
        ------------------------------ \\
        RETURN: \\
        ---------------  \\
        atoms: list -> List Of Element Of All Atoms \\
        ------------------------------ \\
        '''
        # XX AJ please ignore this function for now, I think I'm only using that when I get to the singlepoints on the other branch
        atoms = []
        with open('mass_map.json', 'r') as file:
            mass_map = json.load(file)
          
        name_map = {value: key for key, value in mass_map.items()}
        # XX AJ I need to change qmmm_topology in new version
        qmmm_topology=0
        with open(qmmm_topology) as qmmm_topology_file:
            for line in qmmm_topology_file:
                match = re.search(r"\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    # logger(logfile, "moleculetype section was identified\n")
                    break
            for line in qmmm_topology_file:
                match = re.search(r"\[\s+atoms\s*\]", line, flags=re.MULTILINE)
                if match:
                    # logger(logfile, "atoms section was identified\n")
                    break
            for line in qmmm_topology_file:
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
                        # XX change to logger
                        # if massdiff > 0.01:
                        #     logger(
                        #         logfile,
                        #         str(
                        #             "Found a mass of "
                        #             + str(atommass)
                        #             + " for atom type "
                        #             + str(atomtype)
                        #             + ' (identified as atom name "'
                        #             + str(foundname)
                        #             + '"), which is more than 0.01 different from the expected mass of '
                        #             + str(testmass)
                        #             + ". Atom index was "
                        #             + str(match.group(1))
                        #             + ". This has no effect except unless the atom was identified wrongly or dynamics are intended. Clean your ffnonbonded.itp to avoid these messages!\n"
                        #         ),
                        #     )
                        atoms.append(foundname)
                    else:
                        # logger(
                        #     logfile,
                        #     str(
                        #         "Atom type "
                        #         + str(atomtype)
                        #         + " could not be translated to a regular atom name. Exiting. Last line:\n"
                        #     ),
                        # )
                        # logger(logfile, line)
                        exit(1)
        return atoms


if __name__ == '__main__':

    pass

