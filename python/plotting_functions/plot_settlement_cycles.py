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

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')


def main():
    cycles = np.array([float(10**i) for i in range(7)])

    plt.figure(0, figsize=(12, 16))
    gs1 = gridspec.GridSpec(16, 2)
    plt.figure(1, figsize=(12, 16))
    gs2 = gridspec.GridSpec(16, 2)
    gs = [gs1, gs2]
    axes = [[[], []], [[], []]]
    for i in range(2):
        for j in range(2):
            for fig in range(2):
                plt.figure(fig)
                ax = plt.subplot(gs[fig][i*6:(i+1)*6, j:j+1])
                plt.ylabel('Settlement [mm]', fontsize=24)
                plt.xlabel('Load Cycles [-]', fontsize=24)
                plt.ylim(bottom=0)
                ax.yaxis.set_label_coords(-0.15, 0.5)
                plt.tight_layout()
                axes[fig][i].append(ax)

    plt.draw()
    plt.pause(0.001)
    plt.ion()
    plt.show()

    colors = ['r', 'b', 'g', 'k']
    symbols = ['x', 'o', 's', 'd']
    frequencies = [5., 10., 20., 40.]

    for j, rail_fixture in enumerate(['slab', 'sleepers']):
        for i, geometry in enumerate(['low', 'high']):
            plt.figure(0)
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            ax = plt.subplot(axes[0][i][j])
            max_y1 = 0
            max_y2 = 0
            for load, line in zip([17.5, 22.5, 30.], [':', '-', '--']):
            # for load, line in zip([22.5], ['-']):
                for f, c, sym in zip(frequencies, colors, symbols):
                    odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't_' + str(int(f)) + 'Hz.odb')
                    settlement = 0*cycles
                    for k, n in enumerate(cycles[1:], 1):
                        step_name = 'cycles_' + str(int(n))
                        try:
                            up = abq.get_data_from_path(odb_filename, path_points, 'UP', 'UP2', output_position='NODAL',
                                                        step_name=step_name)
                            settlement[k] = -up[0]*1000
                        except FileNotFoundError:
                            print("Problem with data for {fixture}, {geometry}, {load} tonnes, "
                                  "{f} Hz".format(fixture=rail_fixture, geometry=geometry, load=load, f=f))
                            settlement[k] = np.nan
                    plt.figure(0)
                    if not np.isnan(settlement[-1]):
                        if 10*np.ceil(settlement[-1]/10) > max_y1:
                            max_y1 = 10*np.ceil(settlement[-1]/10)
                            plt.ylim(0, max_y1)
                            plt.draw()
                            plt.pause(0.01)
                    
                    ax.semilogx(cycles[settlement != np.nan], settlement[settlement != np.nan], line + c + sym,
                                lw=2, ms=12, mew=2)
                    if f < 20.:
                        plt.figure(1)
                        ax2 = plt.subplot(axes[1][i][j])
                        cyc = cycles[3:]
                        sett = settlement[3:]
                        sett -= sett[0]
                        ax2.semilogx(cyc, sett, line + c + sym, lw=2, ms=12, mew=2)
                        
                        if 10 * np.ceil(sett[-1] / 10) > max_y2:
                            max_y2 = 10 * np.ceil(sett[-1] / 10)
                            plt.ylim(0, max_y2)
                            plt.draw()
                            plt.pause(0.01)

    fig_lines = [[], []]
    for fig in range(2):
        plt.figure(fig)
        plt.subplot(axes[fig][0][0])
        plt.text(0.05, 0.85, r'\noindent \textbf{Low Embankment\\Concrete Slab}', transform=axes[fig][0][0].transAxes,
                 ma='left', fontweight='bold')

        plt.subplot(axes[fig][0][1])
        plt.text(0.05, 0.85, r'\noindent \textbf{Low Embankment\\Sleepers}', transform=axes[fig][0][1].transAxes,
                 ma='left', fontweight='bold')

        plt.subplot(axes[fig][1][0])
        plt.text(0.05, 0.85, r'\noindent \textbf{High Embankment\\Concrete Slab}', transform=axes[fig][1][0].transAxes,
                 ma='left', fontweight='bold')

        plt.subplot(axes[fig][1][1])
        plt.text(0.05, 0.85, r'\noindent \textbf{High Embankment\\Sleepers}', transform=axes[fig][1][1].transAxes,
                 ma='left', fontweight='bold')
        lines = [
            plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Axle load}')[0],
            plt.plot([0, -1], [-1, -1], ':k', lw=2, label='17.5 t')[0],
            plt.plot([0, -1], [-1, -1], 'k', lw=2, label='22.5 t')[0],
            # plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0],
            plt.plot([0, -1], [-1, -1], '--k', lw=2, label='30.0 t')[0],
            plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{Frequency}')[0]
        ]
        fig_lines[fig] = lines
    plt.figure(0)
    lines = fig_lines[0]
    for f, c, sym, in zip(frequencies, colors, symbols):
        if f in [20]:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0])
        lines.append(plt.plot([0, -1], [-1, -1], '-' + c + sym, lw=2, label=str(int(f)) + ' Hz',
                              ms=12, mew=2)[0])

    leg = plt.legend(handles=lines, ncol=3, bbox_to_anchor=(-1.2, -0.2), loc='upper left',
                     numpoints=1)
    leg.get_texts()[7].set_color("white")
    # leg.get_texts()[2].set_color("white")
    # leg.get_texts()[6].set_color("white")
    plt.ioff()
    plt.savefig(figure_directory + '/settlement_cycles.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"})

    plt.figure(1)
    lines = fig_lines[1]
    for f, c, sym, in zip(frequencies[:2], colors, symbols):
        if f in [20]:
            lines.append(plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'white')[0])
        lines.append(plt.plot([0, -1], [-1, -1], '-' + c + sym, lw=2, label=str(int(f)) + ' Hz',
                              ms=12, mew=2)[0])

    plt.legend(handles=lines, ncol=2, bbox_to_anchor=(-1.0, -0.2), loc='upper left',  numpoints=1)
    plt.ioff()
    plt.savefig(figure_directory + '/settlement_cycles_precomp.tif', dpi=600, pil_kwargs={"compression": "tiff_lzw"})

    plt.show()


if __name__ == '__main__':
    main()
