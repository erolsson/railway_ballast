from __future__ import print_function
import odbAccess

import os

from abaqus_functions.create_empty_odb import create_empty_odb
from abaqus_functions.odb_io_functions import read_field_from_odb

if __name__ == '__main__':
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
        stress = read_field_from_odb(field_id='S', odb_file_name=stress_odb_filename, instance_name=instance_name)
        print(stress.shape)
