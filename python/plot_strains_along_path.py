import os

import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from plot_stresses_along_path import get_tensor_from_path, mises

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/python/embankment_model')
frequency = 10


def main():
    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                ep = get_tensor_from_path(odb_filename, path_points, 'EP')
                e_eff = mises(ep)
                e_vol = ep[:, 0] + ep[:, 1] + ep[:, 2]
                plt.figure(0)
                plt.plot(path_points[0, 1] - path_points[:, 1], e_eff, c + line, lw=2)
                plt.figure(1)
                plt.plot(path_points[0, 1] - path_points[:, 1], -e_vol, c + line, lw=2)
    plt.show()


if __name__ == '__main__':
    main()
