from __future__ import print_function
from collections import namedtuple
import pickle
from common import numpy as np

from abaqusConstants import NODAL
import odbAccess

import os

from abaqus_functions.odb_io_functions import read_field_from_odb
from abaqus_functions.odb_io_functions import write_field_to_odb

simulation_directory = os.path.expanduser('~/railway_ballast/abaqus2014/')
results_odb_filename = simulation_directory + '/Job-14.odb'

BoundaryCondition = namedtuple('BoundaryCondition', ['set_name', 'type', 'component'])

if __name__ == '__main__':
    boundary_conditions = [BoundaryCondition(set_name='M_SURF_BALLAST2', type='surface', component=2),
                           BoundaryCondition(set_name='SET-CENTER', type='node_set', component=1),
                           BoundaryCondition(set_name='SET-FRONT', type='node_set', component=3),
                           BoundaryCondition(set_name='SET-BACK', type='node_set', component=3)]

    results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
    instance_names = results_odb.rootAssembly.instances.keys()
    instance_names = [name for name in instance_names if 'BALLAST' in name]
    step_names = results_odb.steps.keys()
    step_names = [step_name for step_name in step_names if step_name.lower() not in ['gravity', 'train']]
    results_odb.close()

    for step_name in step_names:
        for instance_name in instance_names:
            permanent_deformation, nodes, _ = read_field_from_odb('U', results_odb_filename, step_name='gravity',
                                                                  instance_name=instance_name,
                                                                  set_name='BALLAST',
                                                                  get_position_numbers=True, position=NODAL)
            permanent_deformation *= 0

            results_odb = odbAccess.openOdb(results_odb_filename, readOnly=True)
            instance = results_odb.rootAssembly.instances[instance_name]
            up = np.zeros(len(instance.nodes)*3)
            bc_dofs = []
            for bc in boundary_conditions:
                set_nodes = []
                if bc.type == 'surface':
                    base = results_odb.rootAssembly.surfaces[bc.set_name]
                elif bc.type == 'node_set':
                    base = results_odb.rootAssembly.nodeSets[bc.set_name]
                idx = base.instances.index(instance)
                nodes = base.nodes[idx]
                print(len(nodes), "in set", bc.set_name)
                for n in nodes:
                    set_nodes.append(3*(n.label - 1) + bc.component - 1)
                    bc_dofs.append(3*(n.label - 1) + bc.component - 1)

            bc_dofs = np.unique(np.array(bc_dofs))

            with open('up.pkl') as pickle_handle:
                up_red = pickle.load(pickle_handle)
            print(np.mean(up_red))
            bc_set = set(bc_dofs)
            j = 0
            for i in range(up.shape[0]):
                if i not in bc_set:
                    up[i] = up_red[j]
                    j += 1

            print(np.mean(up))
            counter = 0
            for i, n in enumerate(nodes):
                permanent_deformation[i, :] = up[3*(n.label-1):3*n.label]
            for i in range(3):
                print(np.min(permanent_deformation[:, i]), np.max(permanent_deformation[:, i]),
                      np.mean(permanent_deformation[:, i]))
            results_odb.close()
            write_field_to_odb(permanent_deformation, 'UP', results_odb_filename, step_name,
                               instance_name=instance_name, position=NODAL, set_name='BALLAST')
