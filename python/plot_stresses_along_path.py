import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from finite_element_model.simulations import simulations
from get_data_from_path import get_data_from_path
from comparison_of_models import get_path_points_for_fem_simulation

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/python/embankment_model')


def get_stress_tensor_from_path(odb_file_name, path_points, step_name=None, frame_number=None):
    stress_components = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
    data = np.zeros((path_points.shape[0], 6))
    for i, component in enumerate(stress_components):
        stress = get_data_from_path(path_points, odb_file_name, 'S', step_name=step_name,
                                    frame_number=frame_number, output_position='INTEGRATION_POINT',
                                    component=component)
        data[:, i] = stress
    return data


def main():
    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['low', 'high']:
            odb_filename = odb_directory + '/stresses_' + rail_fixture + '_' + geometry + '.odb'
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            static_stresses = get_stress_tensor_from_path(odb_filename, path_points, step_name='gravity')
            static_pressure = -np.sum(static_stresses[:, :3], axis=1)/3
            plt.figure(0)
            plt.plot(path_points[0, 1] - path_points[:, 1], static_pressure/1e3, 'k' + line, lw=2)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                if load != 30. and rail_fixture != 'sleepers':
                    step_name = 'cyclic_stresses_' + str(load).replace('.', '_') + 't'
                    s = get_stress_tensor_from_path(odb_filename, path_points, step_name=step_name)
                    von_mises = (np.sum(s[:, :3]**2, axis=1) + 3*np.sum(s[:, 3:]**2, axis=1)
                                 - s[:, 0]*s[:, 1] - s[:, 0]*s[:, 2] - s[:, 1]*s[:, 2])**0.5
                    plt.figure(1)
                    plt.plot(path_points[0, 1] - path_points[:, 1], von_mises/1e3, c + line, lw=2)
        name = rail_fixture[0].upper() + rail_fixture[1:]
        plt.figure(0)
        plt.plot([0, -1], [-1, -1], 'k' + line, lw=2, label=name)

    plt.figure(0)
    ax = plt.gca()
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Static pressure, $p_s$ [kPa]')
    plt.xlim(0, 4.3)
    plt.ylim(0, 30)
    plt.text(0.05, 0.9, '(c)', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('../Figures/pressure_graph.png')
    plt.show()


if __name__ == '__main__':
    main()