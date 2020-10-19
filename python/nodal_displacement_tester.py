from __future__ import division, print_function
import os

from abaqusConstants import NODAL

from permanent_deformations import calculate_nodal_displacements_from_strains
from permanent_deformations import BoundaryCondition

from abaqus_functions.odb_io_functions import read_field_from_odb

import common


def main():
    odb_file_name = os.path.expanduser('~/railway_ballast/FEM_test/test_model.odb')
    boundary_conditions = [BoundaryCondition('X0_NODES', 'node_set', 1),
                           BoundaryCondition('Y0_NODES', 'node_set', 2),
                           BoundaryCondition('Z0_NODES', 'node_set', 3),
                           BoundaryCondition('DISP_NODES', 'node_set', 1, values={6: 0.01})]

    displacement_fem, node_labels, _ = read_field_from_odb('U', odb_file_name, step_name='loading', position=NODAL,
                                                           get_position_numbers=True)
    strains, _, element_labels = read_field_from_odb('E', odb_file_name, step_name='loading', get_position_numbers=True)

    displacement_calc, error = calculate_nodal_displacements_from_strains(odb_file_name, boundary_conditions, 'loading',
                                                                          return_error=True)
    gp = 1
    for i, e in enumerate(element_labels):
        print("element=", e, ", gp=", gp, " strains=", strains[i])
        gp += 1
        if gp == 9:
            gp = 1

    for i in range(len(node_labels)):
        print(node_labels[i], displacement_calc[3*i], displacement_fem[i, 0])
        print(node_labels[i], displacement_calc[3*i+1], displacement_fem[i, 1])
        print(node_labels[i], displacement_calc[3*i+2], displacement_fem[i, 2])
    print(error)


if __name__ == '__main__':
    main()
