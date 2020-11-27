import odbAccess

from abaqusConstants import ELEMENT_NODAL

import pickle
import sys

import numpy as np

from abaqus_functions.odb_io_functions import read_field_from_odb
from utilities import BoundaryCondition
sys.path.append('..')
from FEM_functions.elements import C3D8


if __name__ == '__main__':
    settings_pickle_file = sys.argv[-2]
    results_pickle_file = sys.argv[-1]
    with open(settings_pickle_file, 'rb') as settings_pickle:
        data = pickle.load(settings_pickle)
    par = data['parameter_dict']
    boundary_conditions = data['boundary_conditions']
    strain_odb_file_name = str(par['strain_odb_file_name'])
    instance_name = str(par['instance_name'])
    set_name = str(par['set_name'])
    print(set_name)
    step_name = str(par['step_name'])
    strain_field_var = str(par['strain_field_var'])
    strain, _, element_labels = read_field_from_odb(strain_field_var, strain_odb_file_name, step_name=step_name,
                                                    instance_name=instance_name, get_position_numbers=True,
                                                    set_name=set_name,
                                                    frame_number=par['frame_number'])

    _, node_labels, _ = read_field_from_odb(strain_field_var, strain_odb_file_name,
                                            step_name=step_name, instance_name=instance_name, set_name=set_name,
                                            get_position_numbers=True,
                                            frame_number=par['frame_number'], position=ELEMENT_NODAL)

    node_labels = np.array(node_labels, dtype=int)
    node_labels = np.unique(node_labels)
    dofs = np.array([])
    results_odb = odbAccess.openOdb(strain_odb_file_name, readOnly=True)
    if not instance_name:
        instance = results_odb.rootAssembly.instances[results_odb.rootAssembly.instances.keys()[0]]
    else:
        instance = results_odb.rootAssembly.instances[instance_name]
    nodal_displacements = np.zeros(len(node_labels)*3)
    bc_dofs = []
    bc_vals_dict = {}
    for bc in boundary_conditions:
        if bc.type == 'surface':
            base = results_odb.rootAssembly.surfaces[str(bc.set_name)]
        elif bc.type == 'node_set':
            base = instance.nodeSets[str(bc.set_name)]
        else:
            raise ValueError("type attribute of BoundaryCondition must either be \"surface\" or \"node_set\"")
        nodes = base.nodes
        for n in nodes:
            bc_dofs.append(3*(n.label - 1) + bc.component - 1)
        if bc.values is not None:
            for node_label, value in bc.values.items():
                dof = 3*(node_label - 1) + bc.component - 1
                bc_vals_dict[dof] = value

    bc_dofs = np.unique(np.array(bc_dofs))
    elements = {}
    for e_label in np.unique(element_labels):
        element = instance.elements[e_label - 1]
        element_nodes = [instance.nodes[n - 1] for n in element.connectivity]
        elements[e_label] = C3D8(element_nodes)

    with open(results_pickle_file, 'wb') as results_pickle:
        pickle.dump({'strain': strain, 'elements': elements, 'bc_dofs': bc_dofs, 'bc_vals_dict': bc_vals_dict,
                     'nodal_displacements': nodal_displacements, 'element_labels': element_labels},
                    results_pickle)
