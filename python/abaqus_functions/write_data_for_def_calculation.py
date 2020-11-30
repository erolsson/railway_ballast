from __future__ import print_function, division
import odbAccess

from abaqusConstants import ELEMENT_NODAL

from collections import namedtuple
import pickle
import sys

import numpy as np

from odb_io_functions import read_field_from_odb
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
    step_name = str(par['step_name'])
    strain_field_var = str(par['strain_field_id'])

    _, _, element_labels = read_field_from_odb(strain_field_var, strain_odb_file_name, step_name=step_name,
                                               instance_name=instance_name, get_position_numbers=True,
                                               set_name=set_name,
                                               frame_number=-1)

    _, node_labels, _ = read_field_from_odb(strain_field_var, strain_odb_file_name,
                                            step_name=step_name, instance_name=instance_name, set_name=set_name,
                                            get_position_numbers=True,
                                            frame_number=-1, position=ELEMENT_NODAL)

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
            idx = np.where(node_labels == n.label)[0]
            if len(idx):
                bc_dofs.append(3*idx[0] + bc.component - 1)
        if bc.values is not None:
            for node_label, value in bc.values.items():
                idx = np.where(node_labels == node_label)[0]
                if idx:
                    dof = 3*idx[0] + bc.component - 1
                    bc_vals_dict[dof] = value

    bc_dofs = np.unique(np.array(bc_dofs))
    elements = []
    instance_elements = {}
    instance_nodes = {}
    for e in instance.elements:
        instance_elements[e.label] = e
    for n in instance.nodes:
        instance_nodes[n.label] = n

    Node = namedtuple('Node', ['label', 'coordinates'])
    for e_label in np.unique(element_labels):
        element = instance_elements[e_label]

        abaqus_element_nodes = [instance_nodes[n] for n in element.connectivity]
        element_nodes = []
        for abaqus_node in abaqus_element_nodes:
            idx = np.where(node_labels == abaqus_node.label)[0][0]
            element_nodes.append(Node(label=idx, coordinates=abaqus_node.coordinates))
        elements.append(C3D8(element_nodes))

    with open(results_pickle_file, 'wb') as results_pickle:
        pickle.dump({'elements': elements, 'bc_dofs': bc_dofs, 'bc_vals_dict': bc_vals_dict,
                     'nodal_displacements': nodal_displacements, 'node_labels': node_labels},
                    results_pickle)
