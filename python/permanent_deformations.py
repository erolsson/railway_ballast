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


if __name__ == '__main__':
    results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
    instance_names = results_odb.rootAssembly.instances.keys()
    instance_names = [name for name in instance_names if 'BALLAST' in name]
    step_names = results_odb.steps.keys()
    results_odb.close()

    for step_name in step_names:
        for instance_name in instance_names:
            elements = {}
            permanent_strain, _, element_labels = read_field_from_odb('EP', results_odb_filename, step_name=step_name,
                                                                      instance_name=instance_name,
                                                                      get_position_numbers=True)

            results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
            instance = results_odb.rootAssembly.instances[instance_name]
            for e_label in np.unique(element_labels):
                element = instance.elements[e_label - 1]
                element_nodes = [instance.nodes[n-1] for n in element.connectivity]
                elements[e_label] = C3D8(element_nodes)
            row = np.linspace(0, permanent_strain.shape[0]*6*24)
            col = np.zeros(permanent_strain.shape[0]*6*24)
            values = np.zeros(permanent_strain.shape[0]*6*24)
            print(permanent_strain.shape, len(element_labels))
            strain_line = 0
            for i in range(permanent_strain.shape[0]/8):
                element = elements[element_labels[8*i]]
                for j, gp in enumerate(C3D8.gauss_points):
                    B = element.B(*gp)
                    for comp in range(6):
                        try:
                            dispalcement_comp = np.zeros(24)
                            for k, n in enumerate(element.node_labels):
                                dispalcement_comp[3*k] = 3*(n-1)
                                dispalcement_comp[3*k+1] = 3*(n - 1) + 1
                                dispalcement_comp[3*k + 2] = 3*(n - 1) + 2
                            col[strain_line*24:strain_line*24+24] = dispalcement_comp
                        except ValueError:
                            print(strain_line, len(row), permanent_strain.shape)
                            sdfdsfsdffs
                        strain_line += 1
            print(col)
