#!/usr/bin/env python

#   // INITIAL DESCRIPTION //
"""Short Module Description; Reference To Readme"""

#   // MEATDATA // 
__author__ = 'Florian Anders'
__date__ = '2024-01-09'

#   // IMPORTS //

#   Imports Of Existing Libraries
import datetime
import os

#   Imports From Existing Libraries 


#   // TODOS & NOTES //
#   TODO: This Is An Example To Do
#   NOTE: This Is An Example Note

#   // CLASS & METHOD DEFINITIONS //
class Logger():

    '''
    This Class Handles The Logging Of Important Information From gmx2qmmm
    '''

    #   // INIT //
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
    def log_append(path_file_logging: os.path, str_info_to_log: str) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Appends Information To The Specified Log File (file_log)\\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_info_to_log: str -> Text To Be Written To The Logfile \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        with open(file=path_file_logging, mode='a') as file_log_append:

            file_log_append.write\
                (
                    str('-' * 60) + '\n' + \
                    str(datetime.datetime.now()).split('.')[0] + '\n' \
                    '\n' + \
                    str(str_info_to_log) + '\n' + \
                    str('-' * 60) + '\n\n'
                )
            
            file_log_append.close()


    @staticmethod
    def log_overwrite(path_file_logging, str_info_to_log: str) -> None:

        '''
        ------------------------------ \\
        EFFECT: \\
        --------------- \\
        Writes Information To The Specified Log File; \\
        If The File Exists Already, This Overwrites The Contents \\
        ------------------------------ \\
        INPUT: \\
        --------------- \\
        str_info_to_log: str -> Text To Be Written To The Logfile \\
        ------------------------------ \\
        RETURN: \\
        --------------- \\
        NONE \\
        ------------------------------ \\
        '''

        with open(file=path_file_logging, mode='w') as file_log_overwrite:

            file_log_overwrite.write\
                (
                    str('-' * 60) + '\n' + \
                    str(datetime.datetime.now()).split('.')[0] + '\n' \
                    '\n' + \
                    str(str_info_to_log) + '\n' + \
                    str('-' * 60) + '\n\n'
                )
            
            file_log_overwrite.close()

if __name__ == '__main__':

    #   TODO: Unit Testing Should Follow Below!
    
    pass

    