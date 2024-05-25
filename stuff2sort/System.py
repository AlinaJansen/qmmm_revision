#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
# XX AJ I don't know what to put here because the module only has this one class. What is the difference between the module description and the class description here?
"""Short Module Description; Reference To Readme"""

#   // MEATDATA //
__author__ = 'Alina Jansen'
__date__ = '2024-03-25'

#   // IMPORTS //

#   Imports Of Existing Libraries
import re
import os 
import math
import numpy as np

#   Imports From Existing Libraries

#   Imports Of Custom Libraries
import Generators.GeneratorGeometries as geometry

#   Imports From Custom Libraries
from Logging.Logger import Logger
from Generators import generate_pcf_from_top as make_pcf
from Generators._helper import _flatten, get_linkatoms_ang

#   // TODOS & NOTES //
#   TODO:
#   - Add the option to read atom input from a list (if that's what we want, wait for Florians opinion)
#   - Add logger
#   NOTE:

#   // CLASS & METHOD DEFINITIONS //
class SystemInfo():

    '''
    This Class Reads And Stores Information About The System
    '''

    def __init__(self, input_dict) -> None:

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

        self.input_dict = input_dict

        #   Make a list of all topology files
        self.topology_list = self.getincludelist(self.input_dict['topologyfile'])

        # Read The Different Types Of Molecules In The System
        self.list_of_molecules = self.read_molecules()

        # Read The Charges For All Atoms
        self.charge_vector = []
        for element in self.list_of_molecules:
            self.charge_vector.extend(self.readcharges(element))

        #   Read All Atom Lists
        #   XX AJ currently I'm assuming we're always having separate files for atom indices, only remove this comment when we finally decided so or add the possibility for lists
        self.qmatomslist = self.read_atoms_list(self.input_dict['qmatomslist'])
        self.activeatomslist = self.read_atoms_list(self.input_dict['activeatomslist'])
        if input_dict['useinnerouter']:
            self.inneratomslist = self.read_atoms_list(self.input_dict['inneratomslist'])
            self.outeratomslist = self.read_atoms_list(self.input_dict['outeratomslist'])
        #   If we're not using inner/outer we're keeping the lists empty (currently necessary for some functions but that might change later)
        else:
            self.inneratomslist = []
            self.outeratomslist = []

        #   Read Connectivity
        self.connectivity_list = self.read_connectity_from_topology(self.input_dict['topologyfile'])

        #   Read Initial Geometry
        #   XX AJ I would prefer only one function here independent of the file type and only make that distinction within the function. I'll get back to that later when I'm writing GeneratorGeometries.py
        if self.input_dict['coordinatefile'][-4:] == ".gro":
            #logger(logfile, "Reading geometry (.gro)...\n")
            # XX AJ if we rewrite gro files to g96 files, the following function can be deleted and we can read the g96 file with geometry.readg96 afterwards
            self.initial_geometry = geometry.readgeo(self.input_dict['coordinatefile'])
            # Writing high-precision coordinate file
            # logger(logfile, "Writing high-precision coordinate file...")
            self.write_highprec(self.input_dict['coordinatefile']) # XX temp removed logfile until logfile decision of Florian AJ
            self.input_dict['coordinatefile'] = self.input_dict['jobname'] + ".g96"
            # logger(logfile, "Done.\n")
        elif self.input_dict['coordinatefile'][-4:] == ".g96":
            #logger(logfile, "Reading geometry (.g96)...\n")
            self.initial_geometry = geometry.readg96(self.input_dict['coordinatefile'])

        self.number_of_atoms = int(len(self.initial_geometry)/3)

        #   Create xyzq (Coordinates And Charges) For The Whole System
        #   XX AJ same here, I'll make one function of it
        #   XX AJ technically, I would prefer this xyzq function not in this class, but it's used in the get_linkatoms_ang, so I'll keep it here
        if input_dict['useinnerouter']:
            self.xyzq = geometry.make_xyzq_io(self.initial_geometry, self.charge_vector, self.outeratomslist)
        else:
            self.xyzq = geometry.make_xyzq(self.initial_geometry, self.charge_vector)

        #   Read linkatoms and next order atoms
        self.m1list = self.identify_m1()
        self.m2list = self.identify_m2()

        self.linkatoms = get_linkatoms_ang(self.xyzq, self.qmatomslist, self.m1list, self.connectivity_list, [])
        self.linkcorrlist, self.q1list, self.q2list, self.q3list, self.m3list = self.get_linkcorrlist()




    def read_atoms_list(self, atoms_input):

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------
        Reads Atom Indices
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        atoms_input: str -> Index File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        atom_list: list -> List Of Atom Indices \\
        ------------------------------ \\
        '''

        atom_list = []
        with open(atoms_input) as atom_file:
            for line in atom_file:
                if "[" in line or "]" in line:
                    continue
                atom_index_list = re.findall("\d+", line)
                if atom_index_list:
                    atom_list.extend(map(int, atom_index_list))
        atom_list = sorted(np.array(atom_list).astype(int))
        return atom_list

    def read_molecules(self):

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
        molecule_list: list -> List Of Molecules And Their Amount In The System \\
        ------------------------------ \\
        '''

        molecule_list = []
        with open(self.input_dict['topologyfile']) as ifile:
            found = False
            for line in ifile:
                match = re.search(r"^\[ molecules \]", line, flags=re.MULTILINE)
                if match:
                    found = True
                    break
            if not found:
                #   XX AJ turn into Logging
                print('No "molecules" entry in ' + str(self.top) + " found. Exiting.")
                exit(1)
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                else:
                    match = re.search(r"^(\S+)\s+(\d+)", line, flags=re.MULTILINE)
                    if match:
                        molecule_list.append([match.group(1), match.group(2)])
                    else:
                        #   XX AJ turn into Logging
                        print("Found an incomprehensible line in molecules list. Exiting.")
                        print("Last line was:")
                        print(line)
                        exit(1)
        return molecule_list

    def readcharges(self, molecule_entry):
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads Charges For All Atoms In A Molecule Type \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molvecentry: list -> List With Molecule Name And Amount Of This Molecule \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        charge_vector: list -> List Of Charges For Every Atom In This Molecule Types \\
        ------------------------------ \\
        '''

        charge_vector = []
        current_topology = self.input_dict['topologyfile']
        molecule_name = molecule_entry[0]
        molecule_count = int(molecule_entry[1])
        found = self.checkformol(molecule_name, self.input_dict['topologyfile'])

        if not found:
            for element in self.topology_list:
                found = self.checkformol(molecule_name, element)
                if found:
                    current_topology = element
                    break
        if not found:
            print("No charges found for " + str(molecule_name) + ". Exiting.")
            exit(1)
        with open(current_topology) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        matchstring = r"^\s*" + re.escape(molecule_name)
                        match = re.search(matchstring, line, flags=re.MULTILINE)
                        if match:
                            found = True

                            for line in ifile:
                                match = re.search(
                                    r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                )
                                if match:
                                    break
                            break
                        else:
                            found = False
                            break
                    if found:
                        break
            for line in ifile:
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
                    charge_vector.append(float(match.group(1)))

        charge_vector *= molecule_count

        return charge_vector

    def checkformol(self, molecule_name, topology_file):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Checks If The Molecule (molecule_name) Is Listed In The Topology File \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        correct: bool \\
        ------------------------------ \\
        '''

        with open(topology_file) as ifile:
            correct = False
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s*moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        if line.strip() == '' or line.startswith(';'):
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                correct = True
                            break
                if correct:
                    break

        return correct

    def read_connectity_from_topology(self, topology_file):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        count = 0
        connectivity_list = []
        for element in self.list_of_molecules:
            current_topology = self.get_current_topology(element[0], topology_file)
            for i in range(0, int(element[1])):
                mollength = self.get_mollength_direct(element[0], current_topology)
                connset = self.get_connlist(count, element[0], current_topology)
                connectivity_list += connset
                count += int(mollength)

        return connectivity_list

    def get_current_topology(self, molecule_name, topology_file):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of All Molecules In The System \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        current_topology = topology_file
        found = make_pcf.checkformol(molecule_name, topology_file)

        if not found:
            for element in self.topology_list:
                found = make_pcf.checkformol(molecule_name, element)
                if found:
                    current_topology = element
                    break
        if not found:
            # XX AJ replace with logger
            print("No charges found for " + str(molecule_name) + ". Exiting.")
            exit(1)

        return current_topology

    def get_mollength_direct(self, molecule_name, topology_file):
        
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Amount Of Atoms Of The Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of All Molecules \\
        ------------------------------ \\
        '''

        mollength = 0
        with open(topology_file) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[\s+moleculetype\s*\]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in ifile:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[\s*atoms\s*\]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in ifile:
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
                                                        mollength = int(match.group(1))
                                                    match = re.search(
                                                        r"^\s*\n", line, flags=re.MULTILINE
                                                    )
                                                    if match:
                                                        break
                                            break
                                break
                    break
        return mollength

    def get_connlist(self, offset, molecule_name, topology_file):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Function Reads The Connectivity Of One Molecule (molecule_name) \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        offset: XX ?
        molecule_name: str -> Name Of The Molecule \\
        topology_file: str -> Name Of The Topology File \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        connectivity_list: list -> List Of The Connectivity Of One Molecule \\
        ------------------------------ \\
        '''
                
        connectivity_list = []
        with open(topology_file) as ifile:
            for line in ifile:
                match = re.search(r"^;", line, flags=re.MULTILINE)
                if match:
                    continue
                match = re.search(r"^\[ moleculetype \]", line, flags=re.MULTILINE)
                if match:
                    for line in ifile:
                        match = re.search(r"^;", line, flags=re.MULTILINE)
                        if match:
                            continue
                        else:
                            matchstring = r"^\s*" + re.escape(molecule_name)
                            match = re.search(matchstring, line, flags=re.MULTILINE)
                            if match:
                                for line in ifile:
                                    match = re.search(r"^;", line, flags=re.MULTILINE)
                                    if match:
                                        continue
                                    else:
                                        match = re.search(
                                            r"^\[ bonds \]", line, flags=re.MULTILINE
                                        )
                                        if match:
                                            for line in ifile:
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
                                                            offset
                                                        )
                                                        b = int(match.group(2)) + int(
                                                            offset
                                                        )
                                                        if a > b:
                                                            c = b
                                                            b = a
                                                            a = c
                                                        found = False
                                                        for element in connectivity_list:
                                                            if int(element[0]) == a:
                                                                if int(b) not in np.array(
                                                                    element
                                                                ).astype(int):
                                                                    element.append(int(b))
                                                                    found = True
                                                        if not found:
                                                            connectivity_list.append(
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
        return connectivity_list

    def write_highprec(self, gro_file):

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

        # not tested yet, bc I didn't use a gro file AJ
        filename=str(self.input_dict['jobname'] + ".g96")
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

    def identify_m1(self):

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
        m1list: list -> List Of M1 Atom Indices \\
        ------------------------------ \\
        '''

        m1list = []
        for element in self.qmatomslist:
            bondlist = self.get_bondpartners(element)
            for entry in bondlist:
                if (int(entry) not in np.array(self.qmatomslist).astype(int)) and (
                    int(entry) not in np.array(m1list).astype(int)
                ):
                    m1list.append(entry)
        return m1list

    def get_bondpartners(self, target):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All M1 Atoms For A Specific QM Atom \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        target: int -> Index Of Specific QM Atom \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        partnerlist: list -> List Of M1 Atoms For This Specific QM Atom \\
        ------------------------------ \\
        '''

        partnerlist = []
        for entry in self.connectivity_list:
            found = False

            for i in range(0, len(entry)):
                if int(entry[i]) == int(target):
                    found = True
                    break
            if found:
                if int(entry[0]) == int(target):
                    for i in range(1, len(entry)):
                        if int(entry[i]) not in np.array(partnerlist).astype(int):
                            partnerlist.append(int(entry[i]))
                else:
                    if int(entry[0]) not in np.array(partnerlist).astype(int):
                        partnerlist.append(int(entry[0]))

        return partnerlist

    def identify_m2(self):

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
        m2list: list -> List Of M2 Atom Indices \\
        ------------------------------ \\
        '''

        m2list = []
        for element in self.m1list:
            m2line = []
            bondlist = self.get_bondpartners(element)
            for entry in bondlist:
                if int(entry) not in np.array(self.qmatomslist).astype(int):
                    m2line.append(entry)
            m2list.append(m2line)
        return m2list

    def get_linkcorrlist(self):

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Reads All Q1, Q2, Q3 And M3 Atom Indices \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        final_linkcorrlist: list -> ? \\
        q1list: list -> List Of Q1 Atom Indices \\
        q2list: list -> List Of Q2 Atom Indices \\
        q3list: list -> List Of Q3 Atom Indices \\
        m3list: list -> List Of M3 Atom Indices \\
        ------------------------------ \\
        '''
                
        linkcorrlist = []
        m3list = []
        m4list = []
        q1list = []
        q2list = []
        q3list = []
        # get q1 and m2
        for element in self.m1list:
            q1line = []
            for entry in self.connectivity_list:
                if int(element) in np.array(entry).astype(int):
                    if int(element) != int(entry[0]):
                        if int(entry[0]) in np.array(self.qmatomslist).astype(int):
                            q1line.append(int(entry[0]))
                    else:
                        for i in range(1, len(entry)):
                            if int(entry[i]) in np.array(self.qmatomslist).astype(int):
                                q1line.append(int(entry[i]))
            q1list.append(q1line)
        # get q2
        q1list = list(_flatten(q1list))
        for element in q1list:
            q2line = []
            for conn in self.connectivity_list:
                if int(element) in np.array(conn).astype(int):
                    if (
                        int(element) != int(conn[0])
                        and (int(conn[0]) in np.array(self.qmatomslist).astype(int))
                        and (int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)]))
                    ):
                        q2line.append(int(conn[0]))
                    elif int(element) == int(conn[0]):
                        for i in range(1, len(conn)):
                            if (int(conn[i]) in np.array(self.qmatomslist).astype(int)) and (
                                int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                            ):
                                q2line.append(int(conn[i]))
            q2list.append(q2line)
        # get q3
        for element in q2list:
            q3lineline = []
            for entry in element:
                q3line = []
                for conn in self.connectivity_list:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) in np.array(self.qmatomslist).astype(int))
                            and int(conn[0]) not in np.array(q1list).astype(int)
                            and int(conn[0]) not in np.array([int(x) for x in _flatten(q1list)])
                        ):
                            q3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if (
                                    int(conn[i]) in np.array(self.qmatomslist).astype(int)
                                    and int(conn[i]) not in np.array(q1list).astype(int)
                                    and int(conn[i]) not in np.array([int(x) for x in _flatten(q1list)])
                                ):
                                    q3line.append(int(conn[i]))
                q3lineline.append(q3line)
            q3list.append(q3lineline)
        # get m3
        for element in self.m2list:
            m3lineline = []
            for entry in element:
                m3line = []
                for conn in self.connectivity_list:
                    if int(entry) in np.array(conn).astype(int):
                        if (
                            int(entry) != int(conn[0])
                            and (int(conn[0]) not in np.array(self.qmatomslist).astype(int))
                            and int(conn[0]) not in np.array(self.m1list).astype(int)
                        ):
                            m3line.append(int(conn[0]))
                        elif int(entry) == int(conn[0]):
                            for i in range(1, len(conn)):
                                if int(conn[i]) not in np.array(self.qmatomslist).astype(
                                    int
                                ) and int(conn[i]) not in np.array(self.m1list).astype(int):
                                    m3line.append(int(conn[i]))
                m3lineline.append(m3line)
            m3list.append(m3lineline)
        # get m4
        for element in m3list:
            m4linelineline = []
            for entry in element:
                m4lineline = []
                for stuff in entry:
                    m4line = []
                    for conn in self.connectivity_list:
                        if int(stuff) in np.array(conn).astype(int):
                            if (
                                int(stuff) != int(conn[0])
                                and (int(conn[0]) not in np.array(self.qmatomslist).astype(int))
                                and int(conn[0]) not in np.array(self.m1list).astype(int)
                            ):
                                found = False
                                for morestuff in self.m2list:
                                    if int(conn[0]) in np.array(morestuff).astype(int):
                                        found = True
                                        break
                                if not found:
                                    m4line.append(int(conn[0]))
                            elif int(stuff) == int(conn[0]):
                                for i in range(1, len(conn)):
                                    if int(conn[i]) not in np.array(self.qmatomslist).astype(
                                        int
                                    ) and int(conn[i]) not in np.array(self.m1list).astype(int):
                                        found = False
                                        for morestuff in self.m2list:
                                            if int(conn[i]) in np.array(morestuff).astype(
                                                int
                                            ):
                                                found = True
                                                break
                                        if not found:
                                            m4line.append(int(conn[i]))
                    m4lineline.append(m4line)
                m4linelineline.append(m4lineline)
            m4list.append(m4linelineline)
        # set up link atom charge corr pairs: q1-m2, q1-m3, q2-m2, l-m2, l-m3, l-m4 - l are represented as their m1 counterparts!
        count = 0
        for element in self.m1list:
            linkpairline = []
            for entry in self.m2list[count]:
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(m3list[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            for entry in _flatten(m4list[count]):
                if int(element) < int(entry):
                    linkpairline.append([element, entry])
                else:
                    linkpairline.append([entry, element])
            linkcorrlist.append(linkpairline)
            count += 1
        count = 0
        for element in q1list:
            linkpairline = []
            for stuff in self.m2list[count]:
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            for stuff in list(_flatten(m3list[count])):
                if int(element) < int(stuff):
                    linkpairline.append([element, stuff])
                else:
                    linkpairline.append([stuff, element])
            linkcorrlist.append(linkpairline)
            count += 1
        count = 0
        for element in q2list:
            for entry in element:
                linkpairline = []
                for stuff in self.m2list[count]:
                    if int(entry) < int(stuff):
                        linkpairline.append([entry, stuff])
                    else:
                        linkpairline.append([stuff, entry])
                linkcorrlist.append(linkpairline)
            count += 1
        reshaped_linkcorrlist = np.array(list(_flatten(linkcorrlist))).reshape(-1, 2)
        final_linkcorrlist = []
        for element in reshaped_linkcorrlist:
            found = False
            for i in range(0, len(final_linkcorrlist)):
                if (
                    final_linkcorrlist[i][0] == element[0]
                    and final_linkcorrlist[i][1] == element[1]
                ):
                    found = True
                    break
            if found:
                continue
            final_linkcorrlist.append(element)
        return (
            sorted(final_linkcorrlist, key=lambda l: l[0]),
            q1list,
            q2list,
            q3list,
            m3list,
        )
    
    def getincludelist(self, topology_input):

        '''
        ------------------------------ \\
        EFFECT: \\
        ---------------  \\
        Recursively Looks For Other Topology Files Referred To In The Original File  \\
        And Adds Them To The List Of Topology Files \\
        ------------------------------ \\
        INPUT: \\
        ---------------  \\
        topology_input: str -> Name Of Topology \\
        ------------------------------ \\
        RETURN: \\
        ---------------  \\
        toplist: list -> List Of Topology Files \\
        ------------------------------ \\
        '''
        
        toplist = []
        with open(topology_input) as ifile:
            for line in ifile:
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
                        foundname = os.path.join(self.input_dict['gmxtop_path'], *foundname.strip('/').split('/'))
                        check = os.path.isfile(foundname)
                        if not check:
                            print("File " + foundname + " was not found. Maybe update the gmxpath variable in the script? Exiting.")
                            exit(1)
                    toplist.append(foundname)

                    toplist.extend(self.getincludelist(foundname))
        return toplist
    
    def get_atoms(self, qmmm_topology):
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

