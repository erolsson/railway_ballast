import os

import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from plot_stresses_along_path import get_data_from_path, mises

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
frequency = 5


def main():
    for i, rail_fixture in enumerate(['slab', 'sleepers']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                if load != 30. or rail_fixture != 'sleepers':
                    print('\n========================================================================================')
                    print(rail_fixture, geometry, load)
                    print('========================================================================================')
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                    ep = get_data_from_path(path_points, odb_filename, 'EP', component='EP22',
                                            output_position='INTEGRATION_POINT')
                    plt.figure(i)
                    plt.plot(path_points[0, 1] - path_points[:, 1], ep, '--' + c, lw=2)
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz_commonf.odb')
                    ep = get_data_from_path(path_points, odb_filename, 'EP', component='EP22',
                                            output_position='INTEGRATION_POINT')
                    plt.plot(path_points[0, 1] - path_points[:, 1], ep, ':' + c, lw=2)
    plt.show()


if __name__ == '__main__':
    main()
