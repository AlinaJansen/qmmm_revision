import sys
import os
import numpy as np
from unittest.mock import Mock, patch

gmx2qmmm_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(gmx2qmmm_path)

import gmx2qmmm_execute
from gmx2qmmm.operations.qmmm_job import execute_g16, execute_gmx

# set working directory
wd = 'opt_input'
os.chdir(os.path.join('tests', wd))
machine_precision = 10**-8


class IntegrationTest():

    def __init__(self):
        self.correct_output_path = 'opt_output'

    @staticmethod
    def execute_gmx2qmmm_mocking(directory = None):
        if not directory:
            os.chdir(os.path.join('tests', 'opt_input'))
        else:
            os.chdir(directory)

        execute_g16_mock = Mock()
        execute_gmx_mock = Mock()

        with patch('gmx2qmmm.operations.qmmm_job.execute_g16', execute_g16_mock):
            with patch('gmx2qmmm.operations.qmmm_job.execute_gmx', execute_gmx_mock):
                inputFiles = gmx2qmmm_execute.userInputs()
                gmx2qmmm_execute.gmx2qmmm(inputFiles)
        
        return None
    
    def read_oenergy(self, path=''):
        with open(os.path.join(path, 'oenergy.txt'), 'r') as oenergy:
            lines = oenergy.readlines()
            energies = [i.split('\t') for i in lines[1:]]
            floats = np.array([[float(i.strip('\n')) for i in energies[j]] for j in range(len(energies))])
            values = floats[:, 1:]
            return values
        
    def read_oforces(self, path=''):
        with open(os.path.join(path, 'oforces.txt'), 'r') as oforces:
            lines = oforces.readlines()
            data = [line.split() for line in lines]
            number_of_atoms = data.index(['Step1']) - 2
            data = [line[1:] for line in data]
            data = [row for row in data if row]
            result_array = np.array(data, dtype=float)
            number_of_steps = len(result_array) // number_of_atoms
            columns_per_step = len(result_array[0])
            result_array = result_array.reshape(number_of_steps, number_of_atoms, columns_per_step)
            return result_array

    def compare_oenergy(self):
        calculated_energy = self.read_oenergy()
        correct_energy = self.read_oenergy(path = self.correct_output_path)
        difference = calculated_energy - correct_energy
        threshold = np.linalg.norm(np.ones_like(difference)*machine_precision)
        assert np.linalg.norm(difference) < threshold

    def compare_oforces(self):
        calculated_forces = self.read_oforces()
        correct_forces = self.read_oforces(path = self.correct_output_path)
        difference = calculated_forces - correct_forces
        threshold = np.linalg.norm(np.ones_like(difference)*machine_precision)
        assert np.linalg.norm(difference) < threshold

    def analyze(self):
        self.compare_oenergy()
        self.compare_oforces()
    


if __name__ == '__main__':
    Testing = IntegrationTest()
    Testing.execute_gmx2qmmm_mocking()
    Testing.analyze()