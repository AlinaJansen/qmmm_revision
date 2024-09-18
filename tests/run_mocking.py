import sys
import os
import cProfile
import pstats

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..\..'))
integration_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)
pass
sys.path.append(integration_path)

from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx
from IntegrationTest import IntegrationTest

mocking = True

'''
-----------------------------------------------------------------
Insert path to folder here:
-----------------------------------------------------------------
'''
directory = os.path.join(os.path.dirname(__file__), 'opt_input')
prof_file = 'profiling_bfgs_cp24.prof'
txt_file = 'profiling_bfgs_cp24.txt'
'''
-----------------------------------------------------------------
'''
def run_mocking():
    IntegrationTest.execute_gmx2qmmm_mocking(directory)

def runtime_analysis():
    cProfile.run("IntegrationTest.execute_gmx2qmmm_mocking(directory)", filename=os.path.join(directory, prof_file))
    with open(os.path.join(directory, txt_file), "w") as output_file:
        stats = pstats.Stats(os.path.join(directory, prof_file), stream=output_file)
        stats.sort_stats('cumtime')
        stats.print_stats(100)

    print(f"Profiling analysis saved to {txt_file}")

if mocking:
    run_mocking()
else:
    runtime_analysis()
