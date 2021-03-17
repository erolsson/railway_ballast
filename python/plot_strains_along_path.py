import os

import numpy as np

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

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
frequency = 10


def main():
    for fig in [0, 1]:
        plt.figure(fig)
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label='Model')
    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                if load != 30. or rail_fixture != 'sleepers':
                    print('\n========================================================================================')
                    print(rail_fixture, geometry, load)
                    print('========================================================================================')
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                    ep = get_tensor_from_path(odb_filename, path_points, 'EP')
                    ep_eff = mises(ep)
                    ep_vol = -np.sum(ep[:, :3], axis=1)
                    plt.figure(0)
                    plt.plot(path_points[0, 1] - path_points[:, 1], ep_vol, line + c, lw=2)
                    plt.figure(1)
                    plt.plot(path_points[0, 1] - path_points[:, 1], ep_eff, line + c, lw=2)

        name = rail_fixture[0].upper() + rail_fixture[1:]
        for fig in [0, 1]:
            plt.figure(fig)
            plt.plot([0, -1], [-1, -1], 'k' + line, lw=2, label=name)
    plt.figure(0)
    plt.ylim(-0.01, 0.01)
    plt.ylabel('Volumetric strain')
    plt.figure(1)
    plt.ylim(0, 0.06)
    plt.ylabel('Deviatoric strain')

    for fig in [0, 1]:
        plt.figure(fig)
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'$\quad$')
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label='Axle load')
        plt.plot([0, -1], [-1, -1], 'r', lw=2, label='22.5 t')
        plt.plot([0, -1], [-1, -1], 'b', lw=2, label='30 t')
        plt.xlim(0, 4.3)
        plt.xlabel('Distance from ballast surface [m]', fontsize=24)
        plt.legend(loc='best', framealpha=0.5)

    plt.figure(0)
    plt.tight_layout()
    plt.savefig('../Figures/volumetric_strain_graph.png')

    plt.figure(1)
    plt.tight_layout()
    plt.savefig('../Figures/deviatoric_strain_graph.png')
    plt.show()


if __name__ == '__main__':
    main()
