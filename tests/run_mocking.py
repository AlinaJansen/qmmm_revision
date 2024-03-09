import sys
import os

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..'))
integration_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)
sys.path.append(integration_path)

# from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx


'''
-----------------------------------------------------------------
Insert path to folder here:
-----------------------------------------------------------------
'''
directory = os.path.join(os.path.dirname(__file__), 'temp_bfgs_cp24')
'''
-----------------------------------------------------------------
'''

from IntegrationTest import IntegrationTest

IntegrationTest.execute_gmx2qmmm_mocking(directory)
pass


# Next step: untar step2!