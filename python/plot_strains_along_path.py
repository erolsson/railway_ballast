import os

import numpy as np

import matplotlib.gridspec as gridspec
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
    plt.figure(0, figsize=(12, 8))
    gs = gridspec.GridSpec(3, 2)
    ax1 = plt.subplot(gs[0:2, 0:1])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Volumetric strain [-]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(-0.005, 0.01)
    ax1.yaxis.set_label_coords(-0.15, 0.5)
    plt.tight_layout()

    ax2 = plt.subplot(gs[0:2, 1:2])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Deviatoric strain [-]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(0, 0.08)
    ax2.yaxis.set_label_coords(-0.15, 0.5)
    plt.tight_layout()

    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['low', 'high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([22.5, 30.], ['r', 'b']):
                print('\n========================================================================================')
                print(rail_fixture, geometry, load)
                print('========================================================================================')
                odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                ep = get_tensor_from_path(odb_filename, path_points, 'EP')
                ep_eff = mises(ep)
                ep_vol = -np.sum(ep[:, :3], axis=1)
                ax1.plot(path_points[0, 1] - path_points[:, 1], ep_vol, line + c, lw=2)
                ax2.plot(path_points[0, 1] - path_points[:, 1], ep_eff, line + c, lw=2)

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Model}')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='Sleepers')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='Slab')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], 'r', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], 'b', lw=2, label='30 t')[0]
    ]
    _ = plt.legend(handles=lines, ncol=2, bbox_to_anchor=(-0.8, -0.2), loc='upper left')
    plt.savefig('../Figures/strain_graphs.png')
    plt.show()


if __name__ == '__main__':
    main()
