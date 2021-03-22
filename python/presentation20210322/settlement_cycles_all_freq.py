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


def main():
    cycles = np.array([float(10**i) for i in range(3, 7)])

    plt.figure(0, figsize=(12, 8))
    gs = gridspec.GridSpec(3, 1)
    axes = []
    for j in range(1):
        ax = plt.subplot(gs[0:2, 0])
        plt.ylabel('Settlement [mm]', fontsize=24)
        plt.xlabel('Load Cycles [-]', fontsize=24)
        plt.ylim(bottom=0)
        ax.yaxis.set_label_coords(-0.05, 0.5)
        plt.tight_layout()
        axes.append(ax)
        plt.xlim(1e3, 1e6)
    plt.draw()
    plt.pause(0.001)
    plt.ion()
    plt.show()

    colors = ['r', 'b', 'g', 'k']
    symbols = ['x', 'o', 's', 'd']
    frequencies = [5., 10., 20., 40.]

    for j, rail_fixture in enumerate(['sleepers']):
        for i, geometry in enumerate(['low']):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            ax = plt.subplot(axes[j])
            max_y = 0
            for load, line in zip([22.5, 30.], ['-', '--']):
                for f, c, sym in zip(frequencies, colors, symbols):
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(int(f)) + 'Hz.odb')
                    settlement = 0*cycles
                    for k, n in enumerate(cycles):
                        step_name = 'cycles_' + str(int(n))
                        try:
                            up = get_data_from_path(path_points, odb_filename, 'UP', 'UP2', output_position='NODAL',
                                                    step_name=step_name)
                            settlement[k] = -up[0]*1000
                        except FileNotFoundError as err:
                            print(err)
                            print("Problem with data for {fixture}, {geometry}, {load} tonnes, "
                                  "{f} Hz".format(fixture=rail_fixture, geometry=geometry, load=load, f=f))
                            settlement[k] = np.nan

                    y = settlement[settlement != np.nan] - settlement[0]
                    ax.semilogx(cycles[settlement != np.nan], y, line + c + sym, lw=2, ms=12, mew=2)
                    plt.draw()
                    plt.pause(0.01)
                    if not np.isnan(y[-1]):
                        if 10*np.ceil(y[-1]/10) > max_y:
                            max_y = 10*np.ceil(y[-1]/10)
                            plt.ylim(0, max_y)
                            plt.draw()
                            plt.pause(0.01)

    plt.subplot(axes[0])
    plt.text(0.05, 0.85, r'\noindent \textbf{Low Embankment\\Sleepers}', transform=axes[0].transAxes,
             ma='left', fontweight='bold')

    lines = [
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
        plt.plot([0, -1], [-1, -1], 'k', lw=2, label='22.5 t')[0],
        plt.plot([0, -1], [-1, -1], '--k', lw=2, label='30 t')[0],
        plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Frequency}')[0]
    ]
    for f, c, sym, in zip(frequencies, colors, symbols):
        if f in [20]:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0])
        lines.append(plt.plot([0, -1], [-1, -1], '-' + c + sym, lw=2, label=str(int(f)) + ' Hz',
                              ms=12, mew=2)[0])

    leg = plt.legend(handles=lines, ncol=3, bbox_to_anchor=(0.08, -0.2), loc='upper left',
                     numpoints=1)
    leg.get_texts()[6].set_color("white")
    plt.savefig('settlement_1000_cycles_all_freq.png')
    plt.ioff()
    plt.show()


if __name__ == '__main__':
    main()
