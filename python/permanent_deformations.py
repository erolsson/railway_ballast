from __future__ import print_function
from collections import defaultdict
from common import scipy
import pickle
from common import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr

import odbAccess
from abaqusConstants import NODAL

import os
import sys

from abaqus_functions.odb_io_functions import read_field_from_odb
from FEM_functions.elements import C3D8

simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/Job-14.odb'


class BoundaryCondition:
    def __init__(self, set_name, bc_type, component, values=None):
        self.set_name = set_name
        self.type = bc_type
        self.component = component
        self.values = values


def calculate_nodal_displacements_from_strains(odb_file_name, boundary_conditions, step_name, instance_name=None,
                                               strain_field_var='E', frame_number=-1, return_error=False):
    elements = {}
    strain, _, element_labels = read_field_from_odb(strain_field_var, odb_file_name, step_name=step_name,
                                                    instance_name=instance_name,
                                                    get_position_numbers=True, frame_number=frame_number)
    results_odb = odbAccess.openOdb(odb_file_name, readOnly=True)
    if instance_name is None:
        instance = results_odb.rootAssembly.instances[results_odb.rootAssembly.instances.keys()[0]]
    else:
        instance = results_odb.rootAssembly.instances[instance_name]
    nodal_displacements = np.zeros(len(instance.nodes)*3)
    bc_dofs = []
    bc_vals_dict = {}
    for bc in boundary_conditions:
        set_nodes = []
        if bc.type == 'surface':
            base = results_odb.rootAssembly.surfaces[bc.set_name]
        elif bc.type == 'node_set':
            base = instance.nodeSets[bc.set_name]
        else:
            raise ValueError("type attribute of BoundaryCondition must either be \"surface\" or \"node_set\"")
        nodes = base.nodes
        for n in nodes:
            set_nodes.append(3*(n.label - 1) + bc.component - 1)
            bc_dofs.append(3*(n.label - 1) + bc.component - 1)
        if bc.values is not None:
            for node_label, value in bc.values.items():
                dof = 3*(node_label - 1) + bc.component - 1
                bc_vals_dict[dof] = value

    bc_dofs = np.unique(np.array(bc_dofs))
    bc_vals = 0.*bc_dofs
    for i, dof in enumerate(bc_dofs):
        bc_vals[i] = bc_vals_dict.get(dof, 0.)
    nodal_displacements[bc_dofs] = bc_vals
    for e_label in np.unique(element_labels):
        element = instance.elements[e_label - 1]
        element_nodes = [instance.nodes[n - 1] for n in element.connectivity]
        elements[e_label] = C3D8(element_nodes)
    row = np.zeros(strain.shape[0]*6*24)
    col = np.zeros(strain.shape[0]*6*24)
    values = np.zeros(strain.shape[0]*6*24)
    strain_line = 0
    displacement_comp = np.zeros(24)
    ep = np.zeros(strain.shape[0]*6)
    for i in range(strain.shape[0]/8):
        element = elements[element_labels[8*i]]
        for k, n in enumerate(element.node_labels):
            displacement_comp[3*k] = 3*(n - 1)
            displacement_comp[3*k + 1] = 3*(n - 1) + 1
            displacement_comp[3*k + 2] = 3*(n - 1) + 2
        for j, gp in enumerate(C3D8.gauss_points):
            B = element.B(*gp)
            for comp in range(6):
                ep[i*8*6 + j*6 + comp] = strain[i*8 + j, comp]
                col[strain_line*24:strain_line*24 + 24] = displacement_comp
                row[strain_line*24:strain_line*24 + 24] = strain_line
                values[strain_line*24:strain_line*24 + 24] = B[comp, :]
                strain_line += 1
    # col = col[values != 0]
    # row = row[values != 0]
    # values = values[values != 0]
    B_matrix = coo_matrix((values, (row, col)),
                          shape=(strain.shape[0]*6, nodal_displacements.shape[0])).tocsc()
    u = read_field_from_odb('U', odb_file_name, step_name=step_name,
                            instance_name=instance_name, position=NODAL,
                            frame_number=frame_number)
    u_fem = np.zeros(3*u.shape[0])
    u_fem[::3] = u[:, 0]
    u_fem[1::3] = u[:, 1]
    u_fem[2::3] = u[:, 2]
    fem_strain = 0*strain
    fem_strain_vec = (B_matrix*u_fem)
    for i in range(6):
        fem_strain[:, i] = fem_strain_vec[i::6]
    all_cols = np.arange(nodal_displacements.shape[0])
    bc_cols = np.where(np.in1d(all_cols, bc_dofs))[0]
    cols_to_keep = np.where(np.logical_not(np.in1d(all_cols, bc_dofs)))[0]

    B_red = B_matrix[:, cols_to_keep]
    nodal_displacements[cols_to_keep, :] = lsqr(B_red, ep - B_matrix[:, bc_cols]*bc_vals, show=True)[0]

    if return_error:
        error = ep - B_matrix*nodal_displacements
        return nodal_displacements, error
    else:
        return nodal_displacements


def main():
    boundary_conditions = [BoundaryCondition(set_name='M_SURF_BALLAST2', bc_type='surface', component=2),
                           BoundaryCondition(set_name='SET-CENTER', bc_type='node_set', component=1),
                           BoundaryCondition(set_name='SET-FRONT', bc_type='node_set', component=3),
                           BoundaryCondition(set_name='SET-BACK', bc_type='node_set', component=3)]

    results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
    instance_names = results_odb.rootAssembly.instances.keys()
    instance_names = [name for name in instance_names if 'BALLAST' in name]
    step_names = results_odb.steps.keys()
    step_names = [step_name for step_name in step_names if step_name.lower() not in ['gravity', 'train']]
    results_odb.close()
    for step_name in step_names:
        for instance_name in instance_names:
            up_red = calculate_nodal_displacements_from_strains(results_odb_filename, boundary_conditions, step_name,
                                                                instance_name, 'EP')
            with open('up.pkl', 'w') as pickle_handle:
                pickle.dump(up_red, pickle_handle)

            print(np.max(up_red), np.min(up_red))


if __name__ == '__main__':
    main()










