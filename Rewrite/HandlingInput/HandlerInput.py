#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA // 
__author__ = 'Florian Anders'
__date__ = '2024-01-09'

#   // IMPORTS //

#   Imports Of Existing Libraries
import os
import sys

#   Imports From Existing Libraries 
from collections import defaultdict

#   Imports Of Custom Libraries

#   Imports From Custom Libraries
from Logging.Logger import Logger

#   // TODOS & NOTES //
#   TODO: This Is An Example To Do
#   NOTE: This Is An Example Note

#   // CLASS & METHOD DEFINITIONS //
class FileReader():

    '''
    This Class Reads Specified Files And Provides Methods To Extract Specific Information Written In These Files
    '''

    def __init__(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        NONE \\
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

        pass


    @staticmethod
    def read_file(path_file_logging: str, str_filename_input: str) -> defaultdict:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        NONE \\
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

        try:

            with open(str_filename_input) as file_input:

                list_content_file: list = file_input.readlines()

                file_input.close()

            #   Initialise An Empty Defaultdict;
            #   This Will Hold Key / Value Pairs Based On The Inputs
            defaultdict_parameters_input: defaultdict = defaultdict()

            for str_pair_parameter in list_content_file:

                #   Ignore Newline Characters ('\n');
                #   Ignore Empty Lines;
                #   Ignore Lines Beginning With A Hashtag Symbol (#) From The List Of Parameters;
                #   For Now, Ignore Lines Beginning With Exclamation Mark (!); This Is Used For Include Statements, Which We Handle Later
                if str_pair_parameter == '\n' or str_pair_parameter == '' or str_pair_parameter.startswith('#'):

                    continue

                else:
                    
                    #   Populate The Defaultdict For Each Line In The Parameter Input File With:
                    #   str_pair_parameter.strip().split('=')[0] -> The First Part Of The Line Text, Before The Equality Sign ('=') - This Becomes The Key Of The Defaultdict
                    #   str_pair_parameter.strip().split('=')[1] -> The Second Part Of The Line Text, After The Equality Sign ('=') - This Becomes The Value Of The Defaultdict
                    #   See The Syntax Of Constructing A Defaultdict For Further Information On The Following Line Of Code
                    defaultdict_parameters_input[str_pair_parameter.strip().split('=')[0]] = str_pair_parameter.strip().split('=')[1]

            #   Check For Usage Of 'Library' Files After We Have Populated The Initial Defaultdict;
            #   These Should be Included In The Parameter File With '!INCLUDE';
            #   Use The Specified Parameters In The Include File And Override The Previously Created Defaultdict Values
            for str_pair_parameter in list_content_file:

                if str_pair_parameter.startswith('!INCLUDE'):
                    
                    #   Check, If A Value Was Set
                    if str_pair_parameter.strip().split('=')[1] != '':
                        
                        #   If Value Was Set, it Has To Be The File Path + Name (Absolute Path) Of The File To Be Included
                        file_parameters_include = str_pair_parameter.strip().split('=')[1]

                        with open(os.path.abspath(file_parameters_include)) as file_parameters_superior:

                            list_content_file_include: list = file_parameters_superior.readlines()

                            file_parameters_superior.close()
                        
                        #   Iterate Over The List Of Superioir Parameters
                        for str_pair_parameter_superior in list_content_file_include:
                            
                            #   Check, If For Each parameter A Value Has Been Set
                            if str_pair_parameter_superior.split('=')[1] != '':
                                
                                #   Override The Respective Parameter In The Defaultdict With The Value To be Used
                                defaultdict_parameters_input[str_pair_parameter_superior.strip().split('=')[0]] = str_pair_parameter_superior.strip().split('=')[1]

            #   Check, If The Inputs In The Defaultdict Are:
            #   1) Complete For All Required Inputs,
            #   2) All Required Inputs Have Values Assigned,
            #   3) All Numerical Values Are, In Fact, Numerical (int Should Be int, float Should Be float, etc.)
            #   TODO: Implement Assessment!

            Logger.log_append\
                (
                    path_file_logging=path_file_logging,
                    str_info_to_log='Asserting basic correctness of inputs'
                )
            
            try:

                Asserter.assert_defaultdict_input(defaultdict_parameters_input)

                Logger.log_append\
                (
                    path_file_logging=path_file_logging,
                    str_info_to_log='Basic correctness check COMPLETE'
                )

                return defaultdict_parameters_input

            except AssertionError:

                Logger.log_append\
                (
                    path_file_logging=path_file_logging,
                    str_info_to_log=\
                        'Basic correctness check FAILED;' + '\n'\
                        + 'Please make sure that all necessary parameters are provided in your parameter file!'
                )

                print('An error occured! Check the logfile (\'{0}\') for more information.'.format(path_file_logging))

                sys.exit(1)

        except FileNotFoundError:

            Logger.log_append\
                (
                    path_file_logging=path_file_logging,
                    str_info_to_log=\
                        'An error occured when trying to read the parameter file - \n' + \
                        'Make sure that the file exists and try running gmx2qmmm again.'
                )
            
            #   Console Output
            print\
                (
                    '\n\nAn error occured when trying to read the parameter file - \n' + \
                    'Make sure that the file exists and try running gmx2qmmm again.\n\n'
                )
            
        except ValueError:

            Logger.log_append\
                (
                    path_file_logging=path_file_logging,
                    str_info_to_log=\
                        'An error occured when trying to read the parameter file - \n' + \
                        'We encountered a value mismatch. Please check the correctness of your input types!'
                )


class Asserter:

    '''
    This Class Assesses The Correctness Of Inputs (Currently Only On A Basic Level)
    '''

    def __init__(self) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        NONE \\
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

        pass


    @staticmethod
    def assert_defaultdict_input(defaultdict_to_assert: defaultdict) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Asserts The Basic Correctness Of Inputs Stored in A Defaultdict \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        defaultdict_to_assert: defaultdict -> Defaultdict Holding Key / Value Pairs Based On Inputs \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        #   TODO: Expand Assertion Functionality!
        bool_parameters_required_exisitng = False

        list_parameters_input_required =\
            [
                'gaussianpath',
                'gromacspath',
            ]
#gaussianexepath=
#gromacsexepath=
#gromacscoordinatespath=
#gromacstopologypath=
#jobtype=nma
###pccores=3
###pcmemory=
#qmprogram=gaussian
###qmbasisset=
###qmmethod=
#qmatomslist=[1, 2, 3, 4]
###systemcharge=
###systemmultiplicity=
#activeatomslist=[1, 2, 3, 4, 5, 6, 7, 8, 9]
#useinnerouter=False
##inneratomslist=[]
##outeratomslist=[]
        #   TODO: Assertion HERE
        bool_parameters_required_exisitng = True
        if bool_parameters_required_exisitng:

            pass

        else:

            raise AssertionError