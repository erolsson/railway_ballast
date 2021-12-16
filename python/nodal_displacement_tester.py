from __future__ import division, print_function
import os

from deformation_calculator import DeformationCalculator
from abaqus_functions.utilities import BoundaryCondition
from write_data_to_odb import write_data_to_odb
from read_data_from_odb import read_data_from_odb


def main():
    odb_file_name = os.path.expanduser('~/railway_ballast/FEM_test/test_model_second_order.odb')
    boundary_conditions = [BoundaryCondition('X0_NODES', 'node_set', 1),
                           BoundaryCondition('Y0_NODES', 'node_set', 2),
                           BoundaryCondition('Z0_NODES', 'node_set', 3),
                           BoundaryCondition('DISP_NODES', 'node_set', 1, values={6: 0.01})]

    calculator = DeformationCalculator(odb_file_name, boundary_conditions, step_name='loading')
    u, error = calculator.calculate_deformations(step_name='loading')
    u_true = read_data_from_odb('U', odb_file_name, step_name='loading', position='NODAL')
    for i in range(u.shape[0]):
        for axis, j in zip('xyz', range(3)):
            print('D{dof}{axis}\t{u_fem}\t{u_true}\t{u_error}'.format(dof=i+1, axis=axis, u_fem=u[i, j],
                                                                      u_true=u_true[i, j],
                                                                      u_error=u[i, j] - u_true[i, j]))
    sddsddssddsdsds
    write_data_to_odb(u, 'UP', odb_file_name, step_name='loading', position='NODAL', frame_number=1)
    write_data_to_odb(error, 'ERR', odb_file_name, step_name='loading', frame_number=1)


if __name__ == '__main__':
    main()
