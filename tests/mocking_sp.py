import sys
import os
import shutil
from unittest.mock import Mock, patch

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)

import gmx2qmmm_execute
from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx

# set working directory
wd = 'sp_input'
os.chdir(os.path.join('tests', wd))


class IntegrationTest():

    def __init__(self):
        pass

    def execute_gmx2qmmm_mocking():
        execute_g16_mock = Mock()
        execute_gmx_mock = Mock()

        with patch('gmx2qmmm.operations.qmmm_job.execute_g16', execute_g16_mock):
            with patch('gmx2qmmm.operations.qmmm_job.execute_gmx', execute_gmx_mock):
                inputFiles = gmx2qmmm_execute.userInputs()
                gmx2qmmm_execute.gmx2qmmm(inputFiles)
        
        return None
    
IntegrationTest.execute_gmx2qmmm_mocking()