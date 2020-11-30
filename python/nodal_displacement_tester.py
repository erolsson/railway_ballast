from __future__ import division, print_function
import os

from calculate_permanent_deformations import DeformationCalculator
from abaqus_functions.utilities import BoundaryCondition
from write_data_to_odb import write_data_to_odb
from read_data_from_odb import read_data_from_odb


def main():
    odb_file_name = os.path.expanduser('~/railway_ballast/FEM_test/test_model.odb')
    boundary_conditions = [BoundaryCondition('X0_NODES', 'node_set', 1),
                           BoundaryCondition('Y0_NODES', 'node_set', 2),
                           BoundaryCondition('Z0_NODES', 'node_set', 3),
                           BoundaryCondition('DISP_NODES', 'node_set', 1, values={6: 0.01})]

    calculator = DeformationCalculator(odb_file_name, boundary_conditions, step_name='loading')
    u, error = calculator.calculate_deformations(step_name='loading')
    print(u)
    u_true = read_data_from_odb('U', odb_file_name, step_name='loading', position='NODAL')
    print(u_true)
    print(error)
    sddsddssddsdsds
    write_data_to_odb(u, 'UP', odb_file_name, step_name='loading', position='NODAL', frame_number=1)
    write_data_to_odb(error, 'ERR', odb_file_name, step_name='loading', frame_number=1)


if __name__ == '__main__':
    main()
