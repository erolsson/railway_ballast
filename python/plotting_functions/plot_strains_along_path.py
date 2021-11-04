import os

import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.style

from abaqus_python_interface import ABQInterface
from comparison_of_models import get_path_points_for_fem_simulation
from plot_stresses_along_path import mises

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')
abq = ABQInterface("abq2018")
frequency = 5


def main():
    plt.figure(0, figsize=(12, 9))
    gs = gridspec.GridSpec(3, 2)
    ax1 = plt.subplot(gs[0:2, 0:1])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Volumetric strain [-]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(-0.005, 0.01)
    plt.text(0.07, 0.92, r'\bf{(a)}', transform=ax1.transAxes)
    ax1.yaxis.set_label_coords(-0.15, 0.5)
    plt.tight_layout()

    ax2 = plt.subplot(gs[0:2, 1:2])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Deviatoric strain [-]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(0, 0.08)
    plt.text(0.07, 0.92, r'\bf{(b)}', transform=ax2.transAxes)
    ax2.yaxis.set_label_coords(-0.15, 0.5)
    plt.tight_layout()

    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for geometry in ['high']:
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            for load, c in zip([17.5, 22.5, 30.], ['g', 'r', 'b']):
                print('\n========================================================================================')
                print(rail_fixture, geometry, load)
                print('========================================================================================')
                odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                ep = abq.get_tensor_from_path(odb_filename, path_points, 'EP')
                ep_eff = mises(ep)
                ep_vol = -np.sum(ep[:, :3], axis=1)
                ax1.plot(path_points[0, 1] - path_points[:, 1], ep_vol, line + c, lw=2)
                ax2.plot(path_points[0, 1] - path_points[:, 1], ep_eff, line + c, lw=2)

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Model}')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='Sleepers')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='Slab')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label='White')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], 'g', lw=2, label='17.5 t')[0],
        plt.plot([0, -1], [-1, -1], 'r', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], 'b', lw=2, label='30.0 t')[0]
    ]
    leg = plt.legend(handles=lines, ncol=2, bbox_to_anchor=(-0.8, -0.2), loc='upper left')
    leg.get_texts()[3].set_color("white")
    plt.savefig(figure_directory + 'strain_graphs.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.show()


if __name__ == '__main__':
    main()
