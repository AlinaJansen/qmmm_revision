import sys
import os
import cProfile
import pstats

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)

from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx
from IntegrationTest import IntegrationTest

# set working directory
test_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
prof_file = 'profiling_optimization.prof'
txt_file = 'profiling_optimization.txt'

cProfile.run("IntegrationTest.execute_gmx2qmmm_mocking()", filename=os.path.join(test_path, prof_file))
with open(os.path.join(test_path, txt_file), "w") as output_file:
    stats = pstats.Stats(os.path.join(test_path, prof_file), stream=output_file)
    stats.sort_stats('cumtime')
    stats.print_stats(100)

print(f"Profiling analysis saved to {txt_file}")
