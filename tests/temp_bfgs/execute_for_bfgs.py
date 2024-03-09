import sys
import os
import cProfile
import pstats

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..'))
integration_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)
sys.path.append(integration_path)

from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx
from IntegrationTest import IntegrationTest

current_directory = os.path.dirname(os.path.realpath(__file__))
IntegrationTest.execute_gmx2qmmm_mocking(current_directory)
pass