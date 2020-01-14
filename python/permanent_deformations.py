from __future__ import print_function
from common.scipy.sparse import coo_matrix

import odbAccess

import os

from abaqus_functions.odb_io_functions import read_field_from_odb

simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/results.odb'

if __name__ == '__main__':
    results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
    instance_names = results_odb.rootAssembly.instances.keys()
    step_names = results_odb.steps.keys()
    results_odb.close()
    for step_name in step_names:
        for instance_name in instance_names:
            ep = read_field_from_odb('EP', results_odb_filename, step_name=step_name, instance_name=instance_name)
