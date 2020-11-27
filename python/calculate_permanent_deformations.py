from __future__ import print_function, division
import pickle
import os
import subprocess

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr

from common import abq
from FEM_functions.elements import C3D8
from read_data_from_odb import read_data_from_odb
from write_data_to_odb import write_data_to_odb
from abaqus_functions.utilities import BoundaryCondition

simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/Job-14.odb'


def calculate_nodal_displacements_from_strains(strain_odb_file_name, boundary_conditions,
                                               step_name, instance_name='', set_name='', strain_field_var='E',
                                               frame_number=-1):
    work_directory = os.path.dirname(os.path.abspath(strain_odb_file_name))
    parameter_dict = {'strain_field_var': strain_field_var, 'step_name': step_name, 'instance_name': instance_name,
                      'strain_odb_file_name': strain_odb_file_name, 'frame_number': frame_number, 'set_name': set_name}

    parameter_pickle_file = work_directory + '/parameter_strain_pickle.pkl'
    strain_pickle_file = work_directory + '/strain_pickle.pkl'
    with open(parameter_pickle_file, 'wb') as pickle_file:
        pickle.dump({'parameter_dict': parameter_dict, 'boundary_conditions': boundary_conditions}, pickle_file,
                    protocol=2)
    os.chdir('abaqus_functions')
    job = subprocess.Popen(abq + ' python write_permanent_strain_pickle.py ' + parameter_pickle_file + ' '
                           + strain_pickle_file, shell=True)
    job.wait()
    os.chdir('..')
    with open(strain_pickle_file, 'rb') as pickle_file:
        data = pickle.load(pickle_file, encoding='latin1')
    os.remove(strain_pickle_file)
    os.remove(parameter_pickle_file)

    bc_dofs = data['bc_dofs']
    bc_vals_dict = data['bc_vals_dict']
    nodal_displacements = data['nodal_displacements']
    strain = data['strain']
    element_labels = data['element_labels']
    elements = data['elements']
    bc_vals = 0.*bc_dofs
    for i, dof in enumerate(bc_dofs):
        bc_vals[i] = bc_vals_dict.get(dof, 0.)
    nodal_displacements[bc_dofs] = bc_vals

    row = np.zeros(strain.shape[0]*6*24)
    col = np.zeros(strain.shape[0]*6*24)
    values = np.zeros(strain.shape[0]*6*24)
    strain_line = 0
    displacement_comp = np.zeros(24)
    ep = np.zeros(strain.shape[0]*6)
    for i in range(strain.shape[0]//8):
        element = elements[element_labels[8*i]]
        for k, n in enumerate(element.node_labels):
            displacement_comp[3*k] = 3*(n - 1)
            displacement_comp[3*k + 1] = 3*(n - 1) + 1
            displacement_comp[3*k + 2] = 3*(n - 1) + 2
        for j, gp in enumerate(C3D8.gauss_points):
            B = element.B(*gp)
            for comp in range(6):
                ep[i*8*6 + j*6 + comp] = strain[i*8 + j, comp]
                col[strain_line*24:strain_line*24 + 24] = displacement_comp
                row[strain_line*24:strain_line*24 + 24] = strain_line
                values[strain_line*24:strain_line*24 + 24] = B[comp, :]
                strain_line += 1

    B_matrix = coo_matrix((values, (row, col)),
                          shape=(strain.shape[0]*6, nodal_displacements.shape[0])).tocsc()

    all_cols = np.arange(nodal_displacements.shape[0])
    bc_cols = np.where(np.in1d(all_cols, bc_dofs))[0]
    cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, bc_dofs)))[0]

    B_red = B_matrix[:, cols_to_keep]
    num_cols = B_red.shape[1]
    scale_factors = np.ones(num_cols)
    for j in range(num_cols):
        norm = 0
        for val in B_red[:, j]:
            norm += val[0, 0]**2
        norm = np.sqrt(norm)
        B_red[:, j] /= norm
        scale_factors[j] = norm

    calc_displacements = lsqr(B_red, ep - B_matrix[:, bc_cols]*bc_vals, show=True)[0]
    calc_displacements *= scale_factors
    nodal_displacements[cols_to_keep] = calc_displacements

    error = ep - B_matrix*nodal_displacements
    nodal_displacements = nodal_displacements.reshape((nodal_displacements.shape[0]//3, 3))
    error = error.reshape((error.shape[0]//6, 6))
    return nodal_displacements, error


def main():
    pass


if __name__ == '__main__':
    main()










