import os

import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from get_data_from_path import get_data_from_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')


def main():
    cycles = np.array([float(10**i) for i in range(7)])

    plt.figure(0, figsize=(12, 16))
    gs = gridspec.GridSpec(16, 2)
    axes = [[], []]
    for i in range(2):
        for j in range(2):
            ax = plt.subplot(gs[i*6:(i+1)*6, j:j+1])
            plt.ylabel('Settlement [mm]', fontsize=24)
            plt.xlabel('Load Cycles [-]', fontsize=24)
            plt.ylim(bottom=0)
            ax.yaxis.set_label_coords(-0.15, 0.5)
            plt.tight_layout()
            axes[i].append(ax)
    plt.draw()
    plt.pause(0.001)
    plt.ion()
    plt.show()

    colors = ['r', 'b', 'g', 'k']
    symbols = ['x', 'o', 's', 'd']
    frequencies = [5., 10., 20., 40.]

    for j, rail_fixture in enumerate(['slab', 'sleepers']):
        for i, geometry in enumerate(['low', 'high']):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            ax = plt.subplot(axes[i][j])
            max_y = 0
            for load, line in zip([17.5, 22.5, 30.], [':', '-', '--']):
                for f, c, sym in zip(frequencies, colors, symbols):
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(int(f)) + 'Hz.odb')
                    settlement = 0*cycles
                    for k, n in enumerate(cycles[1:], 1):
                        step_name = 'cycles_' + str(int(n))
                        try:
                            up = get_data_from_path(path_points, odb_filename, 'UP', 'UP2', output_position='NODAL',
                                                    step_name=step_name)
                            settlement[k] = -up[0]*1000
                        except FileNotFoundError:
                            print("Problem with data for {fixture}, {geometry}, {load} tonnes, "
                                  "{f} Hz".format(fixture=rail_fixture, geometry=geometry, load=load, f=f))
                            settlement[k] = np.nan

                    ax.semilogx(cycles[settlement != np.nan], settlement[settlement != np.nan], line + c + sym,
                                lw=2, ms=12, mew=2)
                    plt.draw()
                    plt.pause(0.01)
                    if not np.isnan(settlement[-1]):
                        if 10*np.ceil(settlement[-1]/10) > max_y:
                            max_y = 10*np.ceil(settlement[-1]/10)
                            plt.ylim(0, max_y)
                            plt.draw()
                            plt.pause(0.01)

    plt.subplot(axes[0][0])
    plt.text(0.05, 0.85, r'\noindent \textbf{Low Embankment\\Concrete Slab}', transform=axes[0][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[0][1])
    plt.text(0.05, 0.85, r'\noindent \textbf{Low Embankment\\Sleepers}', transform=axes[0][1].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][0])
    plt.text(0.05, 0.85, r'\noindent \textbf{High Embankment\\Concrete Slab}', transform=axes[1][0].transAxes,
             ma='left', fontweight='bold')

    plt.subplot(axes[1][1])
    plt.text(0.05, 0.85, r'\noindent \textbf{High Embankment\\Sleepers}', transform=axes[1][1].transAxes,
             ma='left', fontweight='bold')

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], ':k', lw=2, label='17.5 t')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='30.0 t')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Frequency}')[0]
    ]
    for f, c, sym, in zip(frequencies, colors, symbols):
        if f in [20]:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0])
        lines.append(plt.plot([0, -1], [-1, -1], '-' + c + sym, lw=2, label=str(int(f)) + ' Hz',
                              ms=12, mew=2)[0])

    leg = plt.legend(handles=lines, ncol=3, bbox_to_anchor=(-1.2, -0.2), loc='upper left',
                     numpoints=1)
    leg.get_texts()[7].set_color("white")
    plt.ioff()
    plt.savefig(figure_directory + '/settlement_cycles.png')

    plt.show()


if __name__ == '__main__':
    main()
