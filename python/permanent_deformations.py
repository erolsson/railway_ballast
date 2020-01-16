from __future__ import print_function
from common import scipy
from common import numpy as np
from scipy.sparse import coo_matrix

import odbAccess
from abaqusConstants import ELEMENT_NODAL

import os

from abaqus_functions.odb_io_functions import read_field_from_odb
from FEM_functions.elements import C3D8

simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/results.odb'

instance_names = ['BALLAST_2', 'BALLAST_1']

if __name__ == '__main__':
    results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
    step_names = results_odb.steps.keys()
    results_odb.close()
    elements = {}
    for step_name in step_names:
        for instance_name in instance_names:
            ep, element_labels, _ = read_field_from_odb('EP', results_odb_filename, step_name=step_name,
                                                        instance_name=instance_name, get_position_numbers=True)
            _, _, node_labels = read_field_from_odb('EP', results_odb_filename, step_name=step_name,
                                                    instance_name=instance_name, position=ELEMENT_NODAL,
                                                    get_position_numbers=True)

            results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
            instance = results_odb.rootAssembly.instances[instance_name]
            for e_label in np.unique(element_labels):
                element = instance.elements[e_label - 1]
                element_nodes = [instance.nodes[n] for n in element.connectivity]
                elements[e_label] = C3D8(element_nodes)
