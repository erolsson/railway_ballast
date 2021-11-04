import os

import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.style

from abaqus_python_interface import ABQInterface

from comparison_of_models import get_path_points_for_fem_simulation

abq = ABQInterface("abq2018")

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')

frequency = 5


def main():
    cycles = [100, 1e4, 1e6]

    plt.figure(0, figsize=(12, 16))
    gs = gridspec.GridSpec(16, 2)
    axes = [[], []]
    for i in range(2):
        for j in range(2):
            ax = plt.subplot(gs[i*6:(i+1)*6, j:j+1])
            plt.xlabel('Distance from ballast surface [m]', fontsize=24)
            plt.ylabel('Vertical displacement [mm]', fontsize=24)
            ax.yaxis.set_label_coords(-0.15, 0.5)
            plt.tight_layout()
            axes[i].append(ax)

    plt.subplot(axes[0][0])
    plt.ylim(-0.005, 10)
    plt.xlim(0, 2.2)
    plt.text(0.35, 0.7, r'\noindent \textbf{Low Embankment\\Concrete Slab}', transform=axes[0][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[0][1])
    plt.ylim(-0.005, 30)
    plt.xlim(0, 2.2)
    plt.text(0.35, 0.7, r'\noindent \textbf{Low Embankment\\Sleepers}', transform=axes[0][1].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][0])
    plt.ylim(-0.005, 10)
    plt.xlim(0, 4.3)
    plt.text(0.35, 0.7, r'\noindent \textbf{High Embankment\\Concrete Slab}', transform=axes[1][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][1])
    plt.ylim(-0.005, 40)
    plt.xlim(0, 4.3)
    plt.text(0.35, 0.7, r'\noindent \textbf{High Embankment\\Sleepers}', transform=axes[1][1].transAxes,
             ma='left', fontweight='bold')

    colors = ['r', 'b', 'g', 'k', 'c', 'm']

    for j, rail_fixture in enumerate(['slab', 'sleepers']):
        for i, geometry in enumerate(['low', 'high']):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            ax = plt.subplot(axes[i][j])
            for load, line in zip([17.5, 22.5, 30.], [':', '-', '--']):
                odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                for n, c in zip(cycles, colors):
                    print('\n====================================================================================')
                    print(rail_fixture, geometry, load, n)
                    print('====================================================================================')
                    step_name = 'cycles_' + str(int(n))
                    try:
                        up = abq.get_data_from_path(odb_filename, path_points, 'UP', 'UP2', output_position='NODAL',
                                                    step_name=step_name)
                        ax.plot(path_points[0, 1] - path_points[:, 1], -up*1000, c + line, lw=2)
                    except FileNotFoundError:
                        pass

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], ':k', lw=2, label='17.5 t')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='30.0 t')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Load Cycles}')[0]
    ]
    for n, c, in zip(cycles, colors):
        lines.append(plt.plot([0, -1], [-1, -1], c, lw=2, label=str(int(n)))[0])

    plt.legend(handles=lines, ncol=2, bbox_to_anchor=(-1.0, -0.2), loc='upper left')

    plt.savefig(figure_directory + 'deformation_graphs.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.show()


if __name__ == '__main__':
    main()
