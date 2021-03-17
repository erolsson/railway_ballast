from __future__ import print_function, division

import os
import pickle
import subprocess

import numpy as np

from material_model.material_model import MaterialModel
from multiprocesser.multiprocesser import multi_processer
from common import abq
from write_data_to_odb import write_data_to_odb
from calculate_permanent_deformations import DeformationCalculator
from abaqus_functions.utilities import BoundaryCondition
from material_model.model_parameters import get_parameters
from common import create_temp_dir_name


def evaluate_permanent_strain_for_gp(material_parameters, cycles, static_stress_state, cyclic_stress_state):
    n = static_stress_state.shape[0]
    permanent_strain = np.zeros((len(cycles), n, static_stress_state.shape[1]))
    for i in range(n):
        model = MaterialModel(material_parameters)
        permanent_strain[:, i, :] = model.update(cycles, cyclic_stress_state[i, :], static_stress_state[i, :])
        if i % 1000 == 0:
            print("Evaluated {i} out of {n} points".format(i=i, n=n))
    permanent_strain[:, :, 3:] *= 2
    return permanent_strain


def calculate_permanent_strains(stress_odb_file_name, strain_odb_file_name, cycles, material_parameters):
    try:
        len(cycles)
    except TypeError:
        cycles = np.array(cycles)
    cycles = np.array(cycles)

    work_directory = create_temp_dir_name(strain_odb_file_name)
    os.makedirs(work_directory)

    if not os.path.isfile(strain_odb_file_name):
        os.chdir('abaqus_functions')
        job = subprocess.Popen(abq + ' python create_empty_odb.py ' + strain_odb_file_name + ' ' + stress_odb_file_name,
                               shell=True)
        job.wait()
        os.chdir('..')

    os.chdir('abaqus_functions')
    static_pickle_file = work_directory + '/static_stresses.pkl'
    cyclic_pickle_file = work_directory + '/cyclic_stresses.pkl'
    job = subprocess.Popen(abq + ' python write_stress_state_pickles.py ' + stress_odb_file_name + ' '
                           + static_pickle_file + ' ' + cyclic_pickle_file,
                           shell=True)
    job.wait()
    os.chdir('..')

    with open(static_pickle_file, 'rb') as static_pickle:
        static_data = pickle.load(static_pickle, encoding='latin1')
        static_stresses = static_data['data']/1e3
        instance_name = str(static_data['instance'])
        element_set_name = str(static_data['element_set'])
    os.remove(static_pickle_file)

    with open(cyclic_pickle_file, 'rb') as cyclic_pickle:
        cyclic_stresses = pickle.load(cyclic_pickle, encoding='latin1')['data']/1e3
    os.remove(cyclic_pickle_file)

    n = static_stresses.shape[0]
    permanent_strains = np.zeros((len(cycles), n, static_stresses.shape[1]))
    num_cpus = 12
    chunksize = n//num_cpus
    indices = [i*chunksize for i in range(num_cpus)]
    indices.append(n)
    job_list = []
    for i in range(num_cpus):
        args_list = [material_parameters, cycles, static_stresses[indices[i]:indices[i+1]],
                     cyclic_stresses[indices[i]:indices[i+1]]]
        job_list.append((evaluate_permanent_strain_for_gp, args_list, {}))
    result = multi_processer(job_list, timeout=7200, cpus=num_cpus)
    for i in range(num_cpus):
        permanent_strains[:, indices[i]:indices[i+1], :] = result[i]
    
    for i, n in enumerate(cycles):
        write_data_to_odb(field_data=permanent_strains[i, :, :], field_id='EP', odb_file_name=strain_odb_file_name,
                          step_name='cycles_' + str(n), instance_name=instance_name, set_name=element_set_name)

    boundary_conditions = [BoundaryCondition('X1_NODES', 'node_set', 1),
                           BoundaryCondition('ballast_bottom_nodes', 'node_set', 2),
                           BoundaryCondition('X_SYM_NODES', 'node_set', 1),
                           BoundaryCondition('Z_SYM_NODES', 'node_set', 3),
                           BoundaryCondition('Z1_NODES', 'node_set', 3)]

    calculator = DeformationCalculator(strain_odb_file_name, boundary_conditions, step_name='cycles_' + str(cycles[0]),
                                       instance_name=instance_name, set_name=element_set_name, strain_field_id='EP')
    for i, n in enumerate(cycles):
        up, err = calculator.calculate_deformations(step_name='cycles_' + str(n))

        write_data_to_odb(up, 'UP', strain_odb_file_name, step_name='cycles_' + str(n), position='NODAL',
                          frame_number=0, set_name='EMBANKMENT_INSTANCE_BALLAST_NODES')
        write_data_to_odb(err, 'ERR', strain_odb_file_name, step_name='cycles_' + str(n), frame_number=0,
                          set_name=element_set_name)
    os.removedirs(work_directory)


def main():
    frequencies = [10., 5., 20., 40.]
    load = 30.
    for f in frequencies:
        sim_name = 'sleepers_low_' + str(load).replace('.', '_') + 't'
        cycles = [1, 10, 100, 1000, 10000, 100000, 1000000]
        stress_odb_filename = os.path.expanduser('~/railway_ballast/odbs/embankment_' + sim_name
                                                 + '.odb')

        strain_odb_filename = os.path.expanduser('~/railway_ballast/odbs/results_' + sim_name
                                                 + '_' + str(int(f)) + 'Hz.odb')
        par = get_parameters(frequency=f)
        calculate_permanent_strains(stress_odb_filename, strain_odb_filename, cycles, par)


if __name__ == '__main__':
    main()
