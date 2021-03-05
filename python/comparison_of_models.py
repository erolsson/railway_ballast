import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from embankment_model.simulations import simulations
from get_data_from_path import get_data_from_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = 'embankment_model/'


def get_path_points_for_fem_simulation(sim_name):
    ballast_start_height = 0
    fem_simulation = simulations[sim_name]
    for layer in fem_simulation.layers:
        if layer.name.lower().startswith('ballast'):
            break
        ballast_start_height += layer.height
    total_height = sum([layer.height for layer in fem_simulation.layers])

    path_points = np.zeros((1000, 3))
    y = np.linspace(total_height, ballast_start_height, 1000)
    path_points[:, 1] = y
    path_points[:, 0] = fem_simulation.track_gauge/2
    path_points[:, 2] += 1e-6

    return path_points


def main():
    frequency = 20
    load = 22.5
    cycles = 1000000

    step_name = 'cycles_' + str(int(cycles))
    for i, geometry in enumerate(['low', 'high']):
        for rail_fixture, c in zip(['slab', 'sleepers'], ['b', 'r']):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            print(path_points)

            for calibration, ltype in zip(['', '_commonf'], ['--', ':']):
                odb_file_name = (odb_directory + 'results_' + rail_fixture + '_' + geometry + '_'
                                 + str(load).replace('.', '_') + 't_' + str(int(frequency)) + 'Hz'
                                 + calibration + '.odb')

                up = get_data_from_path(path_points, odb_file_name, 'UP', 'UP2', output_position='NODAL',
                                        step_name=step_name)
                plt.plot(path_points[1, 0] - path_points[1, :], -up*1000, lw=2)
    plt.show()

if __name__ == '__main__':
    main()
