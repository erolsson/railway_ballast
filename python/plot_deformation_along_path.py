import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from plot_permanent_deformations import get_data_from_path


def main():
    odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
    frequency = 10
    cycles = [10**i for i in range(7)]

    for i, rail_fixture in enumerate(['slab', 'sleepers']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                if load != 30. or rail_fixture != 'sleepers':
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                    for n in cycles:
                        print('\n====================================================================================')
                        print(rail_fixture, geometry, load, n)
                        print('====================================================================================')
                        step_name = 'cycles_' + str(int(n))
                        up = get_data_from_path(path_points, odb_filename, 'UP', 'UP2',
                                                output_position='NODAL',
                                                step_name=step_name)
                        plt.figure(i)
                        plt.plot(path_points[0, 1] - path_points[:, 1], -up*1000, c, lw=2)

    plt.figure(0)
    plt.ylim(0, 5)
    plt.ylabel('Vertical displacement [mm]', fontsize=24)
    plt.figure(1)
    plt.ylim(0, 5)
    plt.ylabel('Vertical displacement [mm]', fontsize=24)

    for fig in [0, 1]:
        plt.figure(fig)
        plt.plot([0, -1], [-1, -1], 'r', lw=2, label='22.5 t')
        plt.plot([0, -1], [-1, -1], 'b', lw=2, label='30 t')
        plt.xlim(0, 4.3)
        plt.xlabel('Distance from ballast surface [m]', fontsize=24)
        plt.legend(loc='best', framealpha=0.5)

    plt.show()


if __name__ == '__main__':
    main()
