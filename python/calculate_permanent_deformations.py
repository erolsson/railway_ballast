from __future__ import print_function, division
import pickle
import os
import subprocess

import numpy as np
from scipy.sparse import coo_matrix
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr, norm

from common import abq
from FEM_functions.elements import C3D8
from read_data_from_odb import read_data_from_odb
from multiprocesser.multiprocesser import multi_processer
simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/Job-14.odb'


class DeformationCalculator:
    def __init__(self, strain_odb_file_name, boundary_conditions, step_name, instance_name='', set_name='',
                 strain_field_id='E'):
        print("Init calculator")
        self.odb_file_name = strain_odb_file_name
        self.work_directory = os.path.dirname(os.path.abspath(self.odb_file_name))
        self.stain_field_id = strain_field_id
        self.instance_name = instance_name
        self.set_name = set_name
        parameter_pickle_file = self.work_directory + '/parameter_strain_pickle.pkl'
        data_pickle_file = self.work_directory + '/strain_pickle.pkl'
        parameter_dict = {'instance_name': self.instance_name, 'strain_odb_file_name': self.odb_file_name,
                          'set_name': self.set_name, 'strain_field_id': self.stain_field_id, 'step_name': step_name}
        with open(parameter_pickle_file, 'wb') as pickle_file:
            pickle.dump({'parameter_dict': parameter_dict, 'boundary_conditions': boundary_conditions}, pickle_file,
                        protocol=2)
        os.chdir('abaqus_functions')
        job = subprocess.Popen(abq + ' python write_data_for_def_calculation.py ' + parameter_pickle_file + ' '
                               + data_pickle_file, shell=True)
        job.wait()
        os.chdir('..')
        with open(data_pickle_file, 'rb') as pickle_file:
            data = pickle.load(pickle_file, encoding='latin1')
        os.remove(data_pickle_file)
        os.remove(parameter_pickle_file)

        print("Handling boundary conditions")
        bc_dofs = data['bc_dofs']
        bc_vals_dict = data['bc_vals_dict']
        self.nodal_displacements = data['nodal_displacements']
        elements = data['elements']
        b_components = data['b_components']
        strain_components = data['strain_components']
        displacement_components = data['displacement_components']
        self.bc_vals = 0.*bc_dofs
        for i, dof in enumerate(bc_dofs):
            self.bc_vals[i] = bc_vals_dict.get(dof, 0.)

        self.nodal_displacements[bc_dofs] = self.bc_vals
        row = np.zeros(b_components)
        col = np.zeros(b_components)
        values = np.zeros(b_components)
        print("Assembling B-matrix")
        job_list = []
        for i, element in enumerate(elements):
            job_list.append((calculate_element_data, [element, i], {}))
        print("Starting multiprocessing of elements")
        b_data = multi_processer(job_list, cpus=12, delay=0., timeout=1e5)
        print("Assembling results")
        for i, data in enumerate(b_data):
            n = data[0].shape[0]
            col[i*n:(i + 1)*n] = data[0]
            row[i*n:(i + 1)*n] = data[1]
            values[i*n:(i + 1)*n] = data[2]

        self.B_matrix = coo_matrix((values, (row, col)),
                                   shape=(strain_components, displacement_components)).tocsc()
        print("Shape of B-matrix:", self.B_matrix.shape)
        all_cols = np.arange(self.nodal_displacements.shape[0])
        self.bc_cols = np.where(np.in1d(all_cols, bc_dofs))[0]
        self.cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, bc_dofs)))[0]
        print("Computing scale factors")
        job_list = []
        self.B_red = self.B_matrix[:, self.cols_to_keep]
        num_cols = self.B_red.shape[1]
        for col in range(num_cols):
            job_list.append((calculate_scale_factor, [self.B_red[:, col]], {}))

        print("Scaling B-matrix")
        self.scale_factors = norm(self.B_red, axis=0)
        scale_array = sp.diags([1./self.scale_factors], offsets=[0])
        self.B_red *= scale_array
        print("Init done")

    def calculate_deformations(self, step_name, frame_number=-1):
        strain = read_data_from_odb(self.stain_field_id, self.odb_file_name, step_name=step_name,
                                    frame_number=frame_number, set_name=self.set_name, instance_name=self.instance_name)
        strain = strain.flatten()

        calc_displacements = lsqr(self.B_red, strain - self.B_matrix[:, self.bc_cols]*self.bc_vals,
                                  show=True)[0]
        calc_displacements /= self.scale_factors
        self.nodal_displacements[self.cols_to_keep] = calc_displacements
        error = strain - self.B_matrix*self.nodal_displacements
        nodal_displacements = self.nodal_displacements.reshape((self.nodal_displacements.shape[0]//3, 3))
        error = error.reshape((error.shape[0]//6, 6))
        return nodal_displacements, error


def calculate_element_data(element, idx):
    displacement_comp = np.zeros(element.dofs)
    n = element.dofs*element.strains_components
    dofs = element.dofs
    col = np.zeros(n)
    row = np.zeros(n)
    values = np.zeros(n)
    strain_line = 0
    for k, n in enumerate(element.node_labels):
        displacement_comp[3*k] = 3*n
        displacement_comp[3*k + 1] = 3*n + 1
        displacement_comp[3*k + 2] = 3*n + 2
    for j, gp in enumerate(element.gauss_points):
        B = element.B(*gp)
        for comp in range(6):
            col[strain_line*dofs:strain_line*dofs + dofs] = displacement_comp
            row[strain_line*dofs:strain_line*dofs + element.dofs] = strain_line + idx*element.strains_components
            values[strain_line*dofs:strain_line*dofs + dofs] = B[comp, :]
            strain_line += 1
    return col, row, values


def calculate_scale_factor(column):
    return norm(column)


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
    elements = data['elements']
    bc_vals = 0.*bc_dofs
    for i, dof in enumerate(bc_dofs):
        bc_vals[i] = bc_vals_dict.get(dof, 0.)
    nodal_displacements[bc_dofs] = bc_vals
    row = np.zeros(strain.shape[0]*6*24)
    col = np.zeros(strain.shape[0]*6*24)
    values = np.zeros(strain.shape[0]*6*24)
    ep = np.zeros(strain.shape[0]*6)
    print("Assembling B-matrix")
    job_list = []
    for i, element in enumerate(elements):
        job_list.append((calculate_element_data, [element, i, strain[i*8:(i+1)*8, :]], {}))
    print("Starting multiprocessing of elements")
    b_data = multi_processer(job_list, cpus=12, delay=0., timeout=1e5)
    print("Assembling results")
    for i, data in enumerate(b_data):
        ep[i*48:(i+1)*48] = data[0]
        col[i*6*24*8:(i+1)*6*24*8] = data[1]
        row[i*6*24*8:(i + 1)*6*24*8] = data[2]
        values[i*6*24*8:(i + 1)*6*24*8] = data[3]

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
    calc_displacements /= scale_factors
    nodal_displacements[cols_to_keep] = calc_displacements

    error = ep - B_matrix*nodal_displacements
    nodal_displacements = nodal_displacements.reshape((nodal_displacements.shape[0]//3, 3))
    error = error.reshape((error.shape[0]//6, 6))
    return nodal_displacements, error


def main():
    pass


if __name__ == '__main__':
    main()










