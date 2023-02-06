from __future__ import print_function, division

import os
import numpy as np

from abaqus_python_interface import ABQInterface
from multiprocesser import multiprocesser

from deformation_calculator import DeformationCalculator
from material_model.model_parameters import get_parameters
from material_model.material_model import MaterialModel

abq = ABQInterface("abq2018", output=True)


class BoundaryCondition(object):
    def __init__(self, set_name, bc_type, component, values=None):
        self.set_name = set_name
        self.type = bc_type
        self.component = component
        self.values = values


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


def calculate_permanent_deformations(stress_odb_file_name, strain_odb_file_name, cycles, material_parameters):
    try:
        len(cycles)
    except TypeError:
        cycles = np.array(cycles)
    cycles = np.array(cycles)

    if not os.path.isfile(strain_odb_file_name):
        print("Creating odb", strain_odb_file_name)
        abq.create_empty_odb_from_odb(strain_odb_file_name, stress_odb_file_name)

    # finding the ballast element set name assuming that the model only contains one instance
    odb_dict = abq.get_odb_as_dict(stress_odb_file_name)
    instance_name, instance_data = next(iter(odb_dict['rootAssembly']['instances'].items()))
    element_sets = instance_data['elementSets']
    element_set_name = None
    for element_set_name in element_sets:
        if 'ballast_elements' in element_set_name.lower():
            break
    print("Performing calculations on the element set", element_set_name)
    print("Reading stress states from", stress_odb_file_name)
    static_stresses = abq.read_data_from_odb('S', stress_odb_file_name, step_name='gravity', set_name=element_set_name,
                                             instance_name=instance_name)/1e3
    loading_stresses = abq.read_data_from_odb('S', stress_odb_file_name, step_name='loading', set_name=element_set_name,
                                              instance_name=instance_name)/1e3
    cyclic_stresses = loading_stresses - static_stresses
    print(np.min(static_stresses[:, 1]))
    print(np.min(loading_stresses[:, 1]))
    print("Evaluating permanent strains")
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
    result = multiprocesser.multi_processer(job_list, timeout=7200, cpus=num_cpus)
    for i in range(num_cpus):
        permanent_strains[:, indices[i]:indices[i+1], :] = result[i]

    print("Writing permanent strains to", strain_odb_file_name)
    for i, n in enumerate(cycles):
        abq.write_data_to_odb(field_data=permanent_strains[i, :, :], field_id='EP', odb_file_name=strain_odb_file_name,
                              step_name='cycles_' + str(n), instance_name=instance_name, set_name=element_set_name)

    boundary_conditions = [BoundaryCondition('X1_NODES', 'node_set', 1),
                           BoundaryCondition('ballast_bottom_nodes', 'node_set', 2),
                           BoundaryCondition('X_SYM_NODES', 'node_set', 1),
                           BoundaryCondition('Z_SYM_NODES', 'node_set', 3),
                           BoundaryCondition('Z1_NODES', 'node_set', 3)]
    print("Evaluating permanent deformations")
    calculator = DeformationCalculator(strain_odb_file_name, boundary_conditions, abq=abq,
                                       step_name='cycles_' + str(cycles[0]), instance_name=instance_name,
                                       set_name=element_set_name, strain_field_id='EP')

    print("Writing permanent deformations to", strain_odb_file_name)
    for i, n in enumerate(cycles):
        up, err = calculator.calculate_deformations(step_name='cycles_' + str(n))

        abq.write_data_to_odb(up, 'UP', strain_odb_file_name, step_name='cycles_' + str(n), position='NODAL',
                              frame_number=0, set_name='EMBANKMENT_INSTANCE_BALLAST_NODES')
        abq.write_data_to_odb(err, 'ERR', strain_odb_file_name, step_name='cycles_' + str(n), frame_number=0,
                              set_name=element_set_name)


def main():
    frequencies = [10., 5., 20., 40]
    load = 17.5
    for f in frequencies:
        sim_name = 'sleepers_low_' + str(load).replace('.', '_') + 't'
        cycles = [1, 10, 100, 1000, 10000, 100000, 1000000]
        stress_odb_filename = os.path.expanduser('~/hassan/embankment_' + sim_name
                                                 + '_2018.odb')

        strain_odb_filename = os.path.expanduser('~/hassan/results_' + sim_name
                                                 + '_' + str(int(f)) + 'Hz.odb')
        par = get_parameters(frequency=f)
        calculate_permanent_deformations(stress_odb_filename, strain_odb_filename, cycles, par)


if __name__ == '__main__':
    main()
