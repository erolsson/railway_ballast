import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from get_data_from_path import get_data_from_path
from comparison_of_models import get_path_points_for_fem_simulation

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')


def get_tensor_from_path(odb_file_name, path_points, field_id, step_name=None, frame_number=None):
    components = ['11', '22', '33', '12', '13', '23']
    data = np.zeros((path_points.shape[0], 6))
    for i, component in enumerate(components):
        stress = get_data_from_path(path_points, odb_file_name, field_id, step_name=step_name,
                                    frame_number=frame_number, output_position='INTEGRATION_POINT',
                                    component=field_id + component)
        data[:, i] = stress
    return data


def mises(tensor):
    return (np.sum(tensor[:, :3]**2, axis=1) + 3*np.sum(tensor[:, 3:]**2, axis=1)
            - tensor[:, 0]*tensor[:, 1] - tensor[:, 0]*tensor[:, 2] - tensor[:, 1]*tensor[:, 2])**0.5


def main():
    plt.figure(1)
    plt.plot([0, -1], [-1, -1], 'w', lw=2, label='Model')
    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                odb_filename = (odb_directory + '/stresses_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't.odb')

                if load == 22.5:
                    static_stresses = get_tensor_from_path(odb_filename, path_points, 'S', step_name='gravity')
                    static_pressure = -np.sum(static_stresses[:, :3], axis=1)/3
                    plt.figure(0)
                    plt.plot(path_points[0, 1] - path_points[:, 1], static_pressure/1e3, 'k' + line, lw=2)

                s = get_tensor_from_path(odb_filename, path_points, 'S', step_name='cyclic_stresses')
                von_mises = mises(s)
                plt.figure(1)
                plt.plot(path_points[0, 1] - path_points[:, 1], von_mises/1e3, c + line, lw=2)
        name = rail_fixture[0].upper() + rail_fixture[1:]
        plt.figure(0)
        plt.plot([0, -1], [-1, -1], 'k' + line, lw=2, label=name)

        plt.figure(1)
        plt.plot([0, -1], [-1, -1], 'k' + line, lw=2, label=name)

    plt.figure(1)
    plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'$\quad$')
    plt.plot([0, -1], [-1, -1], 'w', lw=2, label='Axle load')
    plt.plot([0, -1], [-1, -1], 'r', lw=2, label='22.5 t')
    plt.plot([0, -1], [-1, -1], 'b', lw=2, label='30 t')

    plt.figure(0)
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Static pressure, $p_s$ [kPa]')
    plt.xlim(0, 4.3)
    plt.ylim(0, 30)
    plt.legend(loc='upper left', bbox_to_anchor=[0.6, 0.5])
    plt.tight_layout()
    plt.savefig('../Figures/pressure_graph.png')

    plt.figure(1)
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Cyclic von Mises stress, $q$ [kPa]')
    plt.xlim(0, 4.3)
    plt.ylim(0, 150)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('../Figures/von_mises_graph.png')

    plt.show()


if __name__ == '__main__':
    main()