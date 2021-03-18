import os

import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from plot_permanent_deformations import get_data_from_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
frequency = 10


def main():
    cycles = [10**i for i in range(1, 7)]

    plt.figure(0, figsize=(12, 14))
    gs = gridspec.GridSpec(5, 2)
    axes = [[], []]
    for i in range(2):
        for j in range(2):
            ax = plt.subplot(gs[i*2:(i+1)*2, j:j+1])
            plt.xlabel('Distance from ballast surface [m]', fontsize=24)
            plt.ylabel('Vertical displacement [mm]', fontsize=24)
            ax.yaxis.set_label_coords(-0.15, 0.5)
            plt.tight_layout()
            axes[i].append(ax)

    plt.subplot(axes[0][0])
    plt.ylim(-0.005, 10)
    plt.xlim(0, 2.2)
    plt.text(0.35, 0.7, r'\noindent \textbf{Concrete Slab\\Low Embankment}', transform=axes[0][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[0][1])
    plt.ylim(-0.005, 30)
    plt.xlim(0, 2.2)
    plt.text(0.35, 0.7, r'\noindent \textbf{Sleepers\\Low Embankment}', transform=axes[0][1].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][0])
    plt.ylim(-0.005, 20)
    plt.xlim(0, 4.3)
    plt.text(0.35, 0.7, r'\noindent \textbf{Concrete Slab\\High Embankment}', transform=axes[1][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][1])
    plt.ylim(-0.005, 40)
    plt.xlim(0, 4.3)
    plt.text(0.35, 0.7, r'\noindent \textbf{Sleepers\\High Embankment}', transform=axes[1][1].transAxes,
             ma='left', fontweight='bold')

    colors = ['r', 'b', 'g', 'k', 'c', 'm']

    for j, rail_fixture in enumerate(['slab', 'sleepers']):
        for i, geometry in enumerate(['low', 'high']):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            ax = plt.subplot(axes[i][j])
            for load, line in zip([22.5, 30.], ['-', '--']):
                odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                + str(load).replace('.', '_') + 't_' + str(frequency) + 'Hz.odb')
                for n, c in zip(cycles, colors):
                    print('\n====================================================================================')
                    print(rail_fixture, geometry, load, n)
                    print('====================================================================================')
                    step_name = 'cycles_' + str(int(n))
                    try:
                        up = get_data_from_path(path_points, odb_filename, 'UP', 'UP2', output_position='NODAL',
                                                step_name=step_name)
                        ax.plot(path_points[0, 1] - path_points[:, 1], -up*1000, c + line, lw=2)
                    except FileNotFoundError:
                        pass

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='30 t')[0],
    ]
    for n, c, in zip(cycles, colors):
        if n in [10, 100000]:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0])
        if n == 1000:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Load Cycles}')[0])
        lines.append(plt.plot([0, -1], [-1, -1], c, lw=2, label=str(int(n)))[0])

    leg = plt.legend(handles=lines, ncol=4, bbox_to_anchor=(-1.6, -0.2), loc='upper left')
    leg.get_texts()[3].set_color("white")
    leg.get_texts()[9].set_color("white")
    plt.savefig('../Figures/deformation_graphs.png')
    plt.show()


if __name__ == '__main__':
    main()
