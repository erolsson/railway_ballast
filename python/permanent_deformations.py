from __future__ import print_function
from common import numpy as np
from common import scipy
import odbAccess

import os

from abaqus_functions.create_empty_odb import create_empty_odb
from abaqus_functions.odb_io_functions import read_field_from_odb

from material_model.functional_shapes import permanent_strain
from multiprocesser.multiprocesser import multi_processer

if __name__ == '__main__':
    cycles = 1e5

    material_parameters = np.array([5.19502802e+00, 7.64481869e-07, 4.70091315e+01, 1.34948445e-13,
                                    1.14929617e+01, 1.12439205e+01])

    simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
    stress_odb_filename = simulation_directory + '/embankment_EO.odb'
    results_odb_filename = simulation_directory + '/results.odb'
    if not os.path.exists(results_odb_filename):
        print("Creating new odb with name:")
        print("\t" + results_odb_filename)
        create_empty_odb(results_odb_filename, stress_odb_filename)

    # Looking through the instances to find ballast instances in stress odb
    stress_odb = odbAccess.openOdb(stress_odb_filename, readOnly=True)
    instance_names = stress_odb.rootAssembly.instances.keys()
    for instance_name in instance_names:
        stress = read_field_from_odb(field_id='S', odb_file_name=stress_odb_filename, instance_name=instance_name)/1e6
        pressure = -np.sum(stress[:, 0:3], 1)
        deviator = np.copy(stress)
        for i in range(3):
            deviator[:, i] -= pressure
        von_Mises = np.sqrt(np.sum(stress[:, 0:3]**2, 1) - stress[:, 0]*stress[:, 1] -
                            stress[:, 0]*stress[:, 2] - stress[:, 1]*stress[:, 2] + 3*np.sum(stress[:, 3:]**2, 1))
        direction = 1.5*deviator
        for i in range(6):
            direction[:, i] /= von_Mises
        if "BALLAST" in instance_name:
            print(instance_name)
            job_list = []
            for i in range(pressure.shape[0]):
                job_list.append((permanent_strain, [],
                                 {"cycles": cycles, "p": pressure[i], "q": von_Mises[i],
                                  'parameters': material_parameters}))
            ep = direction*multi_processer(job_list, delay=0., timeout=3600, cpus=8)
            print(np.max(np.max(ep, 1)))
        else:
            ep = 0*deviator
