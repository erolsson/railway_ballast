from __future__ import print_function
from common import numpy as np
import odbAccess

import os

from abaqus_functions.create_empty_odb import create_empty_odb
from abaqus_functions.odb_io_functions import read_field_from_odb
from abaqus_functions.odb_io_functions import write_field_to_odb

from material_model.functional_shapes import permanent_strain
from multiprocesser.multiprocesser import multi_processer


if __name__ == '__main__':
    cycles = 1e5

    material_parameters = np.array([5.19502802e+00, 7.64481869e-07, 4.70091315e+01, 1.34948445e-13,
                                    1.14929617e+01, 1.12439205e+01])

    simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
    stress_odb_filename = simulation_directory + '/Job-14.odb'
    results_odb_filename = simulation_directory + '/results_20200119.odb'
    if not os.path.exists(results_odb_filename):
        print("Creating new odb with name:")
        print("\t" + results_odb_filename)
        create_empty_odb(results_odb_filename, stress_odb_filename)

    # Looking through the instances to find ballast instances in stress odb
    stress_odb = odbAccess.openOdb(stress_odb_filename, readOnly=True)
    instance_names = stress_odb.rootAssembly.instances.keys()
    stress_odb.close()
    for instance_name in instance_names:
        static_stress = read_field_from_odb(field_id='S', odb_file_name=stress_odb_filename,
                                            instance_name=instance_name, step_name='gravity')/1e3
        stress2 = read_field_from_odb(field_id='S', odb_file_name=stress_odb_filename,
                                      instance_name=instance_name, step_name='gravity')/1e3
        ep = 0*static_stress
        if "BALLAST" in instance_name:
            cyclic_stress = stress2 - static_stress
            pressure_static = -np.sum(static_stress[:, 0:3], 1)/3
            pressure_cyclic = -np.sum(cyclic_stress[:, 0:3], 1)/3
            deviator = np.copy(cyclic_stress)
            for i in range(3):
                deviator[:, i] += pressure_cyclic

            von_Mises_cyclic = np.sqrt(np.sum(cyclic_stress[:, 0:3]**2, 1) - cyclic_stress[:, 0]*cyclic_stress[:, 1]
                                       - cyclic_stress[:, 0]*cyclic_stress[:, 2]
                                       - cyclic_stress[:, 1]*cyclic_stress[:, 2]
                                       + 3*np.sum(cyclic_stress[:, 3:]**2, 1))
            print(von_Mises_cyclic[0:10])
            direction = 1.5*deviator
            for i in range(6):
                direction[:, i] /= von_Mises_cyclic

            print(instance_name)
            job_list = []
            for i in range(pressure_static.shape[0]):
                job_list.append((permanent_strain, [],
                                 {"cycles": cycles, "p": pressure_static[i], "q": von_Mises_cyclic[i],
                                 'parameters': material_parameters}))
            ep_magnitude = np.array(multi_processer(job_list, delay=0., timeout=3600, cpus=8))
            ep_magnitude[ep_magnitude > 1.] = 1.
            for i in range(6):
                ep[:, i] = ep_magnitude*direction[:, i]
            max_idx = np.argmax(ep_magnitude)
            print("Maximum permanent strain is", ep[max_idx], "The pressure is", pressure_static[max_idx],
                  "and von Mises is", von_Mises_cyclic[max_idx])

        write_field_to_odb(field_data=ep, field_id='EP', odb_file_name=results_odb_filename,
                           step_name=str(cycles) + '_cycles', instance_name=instance_name, frame_number=0)
