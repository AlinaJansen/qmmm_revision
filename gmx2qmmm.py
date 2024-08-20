#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA // 
__author__ = 'Florian Anders'
__date__ = '2024-01-09'

#   // IMPORTS //

#   Imports Of Existing Libraries
import argparse
import datetime
import os
import textwrap

#   Imports From Existing Libraries
from collections import defaultdict

#   Imports Of Custom Libraries
import System.System as System
import Generators.GeneratorTopologies as Top
import Generators.GeneratorPCF as PCF
import Jobs.Singlepoint as SP

#   Imports From Custom Libraries
from Logging.Logger import Logger
from HandlingInput.HandlerInput import FileReader

#   // TODOS & NOTES //
#   TODO: This Is An Example To Do
#   NOTE: This Is An Example Note
#   TODO: Read And Include GROMACS / VMD Style Index Files; Extractable From Different Programs(?)
#   TODO: Include Prio For Library / Parameter Files

#   // CLASS & METHOD DEFINITIONS //
class GMX2QMMM():

    '''
    This Class Is The Main Object Of Every gmx2qmmm Run
    And Handles The Initial Arguments And Input Files
    '''

    #   // INIT //
    def __init__(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Initialises The gmx2qmmm Members \\
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
        
        #   // DIRECTORIES //
        #   Define The Current Working Directory
        self.directory_base: str = os.path.dirname(os.path.abspath(__file__))

        #   // PARSER //
        #   Declare gmx2qmmm_parser Member For Input File Processing
        self.gmx2qmmm_parser: argparse.ArgumentParser = argparse.ArgumentParser\
            (
                description="gmx2qmmm, a python interface for Quantum mechanics/Molecular mechanics (QM/MM) calculation.",
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=textwrap.dedent\
                    (
                        '''\n
                        gmx2qmmm (GNU licenced) by Jan Philipp Goetze.\n
                        gmx2qmmm is a free QM/MM interface for Gromacs. Compatible with most versions, starting from 5.X.\n
                        gmx2qmmm uses an additive QM/MM scheme with a numerical charge shift procedure and correction of cut bonds via model potentials.\n
                        '''
                    )
            )

        #   Add Arguments (Files) To The Parser (Alphabetical Order Per Subset)
        #   Structure Defining Files

        #   Additional Files
        self.gmx2qmmm_parser.add_argument\
            (
                "-l",
                "--logfile",
                help="logfile(.log)",
                type=str,
                default="logfile.txt"
            )

        self.gmx2qmmm_parser.add_argument\
            (
                "-p",
                "--parameterFile",
                help="Parameter File(.txt)",
                type=str,
                # default="path.dat"
                default="params.txt"    # XX AJ Do you mean params.txt?
            )

        #   Parse And Create Input File Member
        self.files_input_initial: argparse.Namespace = self.gmx2qmmm_parser.parse_args()

        #   Set File Path For Logging (Absolute Path)
        self.path_filepath_logfile = os.path.abspath(self.files_input_initial.logfile)

        #   // LOGGING //
        #   Check, if A Log File Already Exists;
        #   If So, Rename The Old Log File So We Can Freshly Write To A new Logfile Without Losing The Previous One
        if os.path.exists(self.path_filepath_logfile):

            os.rename\
                (
                    self.path_filepath_logfile,
                    '{0}_{1}_{2}.txt'.format\
                        (
                            str(self.path_filepath_logfile).split('.txt')[0],
                            datetime.datetime.now().date(), 
                            str(datetime.datetime.now().time()).split('.')[0].replace(':', '-')
                        )
                )

            Logger.log_append\
                (
                    path_file_logging=self.path_filepath_logfile,
                    str_info_to_log='Initialisation of gmx2qmmm completed; Renamed an older log file!'
                )
        else:

            Logger.log_append\
                (
                    path_file_logging=self.path_filepath_logfile,
                    str_info_to_log='Initialisation of gmx2qmmm completed!'
                )

        #   Create Defaultdict Holding Parameters From input File
        self.defaultdict_parameters_input: defaultdict = FileReader.read_file\
            (
                path_file_logging=self.path_filepath_logfile,
                str_filename_input=self.files_input_initial.parameterFile
            )
        
        #   Format The Defaultdict's Key / Value Pairs So We Can Log Them Nicely;
        #   String Building Time! :D
        str_parameters_input_log: str = ''

        for key, value in self.defaultdict_parameters_input.items():

            str_parameters_input_log += '{0}:{1}{2}\n'.format(key, ' ' * (25 - len(key)), value)

        #   Log The Input Parameters We Just Build
        Logger.log_append\
            (
                path_file_logging=self.path_filepath_logfile,
                str_info_to_log= 'Listing the parameters given as input for this run of gmx2qmmm:\n\n' + str_parameters_input_log
            )



        #   Initializing The System
        self.initalize_system()

        #   Generate QMMM Topology
        self.generate_topology()

        #   Add List Of Elements To System Instance
        self.system.list_atom_elements = self.system.get_atoms(self.topology.qmmm_topology)

        #   Initialize Pointchargefield
        self.generate_PCF()


        #   // JOB SPAWNING //
        #   Assess Job Type
        self.assess_job()

        #   Start Job According To Inputs
        self.start_job()

    def initalize_system(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Initializes The System \\
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
        self.system = System.SystemInfo(self.defaultdict_parameters_input)
        
    def generate_topology(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Writes Large QMMM Topology File \\
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
        self.topology = Top.GenerateTopology(self.defaultdict_parameters_input, self.system, self.directory_base)

    def generate_PCF(self) -> None:
        
        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Writes Large QMMM Topology File \\
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
        # XX AJ I forgot what I used the Job keyword for, I think I will only need it later, I will get back to that
        self.pointchargefield = PCF.GeneratePCF(self.defaultdict_parameters_input, self.system, self.topology, self.directory_base, job=None)

    def assess_job(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        This Function Provides Basic Assessment Of The Type Of Job To Be Performed \\
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

        Logger.log_append\
            (
                path_file_logging=self.path_filepath_logfile,
                str_info_to_log='Assessing job type'
            )

        #   TODO: Implementation Of Assessment Of Basic Correctness!

        Logger.log_append\
            (
                path_file_logging=self.path_filepath_logfile,
                str_info_to_log='Starting requested job: {0}'.format(self.defaultdict_parameters_input['jobtype'])
            )
        
    def start_job(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        This Function Starts The Requested Job \\
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

        if self.defaultdict_parameters_input['jobtype'] == "SINGLEPOINT":
            # logger(logfile, "Performing an single point calculation.\n")
            # SP.Singlepoint(self.defaultdict_parameters_input, self.system, self.topology, self.pcf.pcf_filename)
            SP.Singlepoint(self.defaultdict_parameters_input, self.system, self.topology, self.pointchargefield, self.directory_base)
            pass

        # elif jobtype == "OPT":
        #     # logger(logfile, "Performing an optimization.\n")
        #     # logger(logfile, "Getting initial energy:\n")
        #     perform_opt(qmmmInputs)

        # elif jobtype == "NMA":
        #     jobname = stepper(qmmminfo[0], step)
        #     perform_nma(qmmmInputs)

        # elif jobtype == "SCAN":
        #     # logger(logfile, "Performing scan calculations.\n")
        #     perform_scan(qmmmInputs)

        # elif jobtype == "OPT_ROOTFOLLOWING":
        #     # logger(logfile, "Performing an optimization.\n")
        #     # logger(logfile, "Getting initial energy:\n")
        #     perform_opt_root(qmmmInputs)


#   // MAIN //
if __name__ == "__main__":

    # // INSTANTIATION //
    gmx2qmmm_main = GMX2QMMM()

    #   Print The Defaultdict (Comment Out, Not Needed!)
    #print('defaultdict:\n\n{0}'.format(gmx2qmmm_main.defaultdict_parameters_input))

    #   TODO: Unit Test?