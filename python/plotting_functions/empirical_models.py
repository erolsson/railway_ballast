import os

import numpy as np

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


def selig_2(n, _):
    return 4.318*n**0.17


def shenton_2(n, axle_load):
    ks = 1.1
    return ks*axle_load/20*(0.69*n**0.2 + 2.7e-6*n)


def thom_1(n, _):
    s = (np.log10(n) - 2.4)**2
    s[np.log10(n) < 2.4] = 0
    return s


def hettler(n, axle_load):
    r = 0.00095
    c = 0.43
    f = axle_load*9.82/2
    return r*f**1.6*(1+c*np.log(n))


def main():
    fig = plt.figure(0)
    cycles = np.array([float(10**i) for i in range(7)] + [2e6])
    plt.semilogx(cycles, thom_1(cycles, 0), 'g', lw=2)
    plt.semilogx(cycles, selig_2(cycles, 0), 'k', lw=2)
    path_points = get_path_points_for_fem_simulation('sleepers_high')

    for load, line in zip([17.5, 22.5, 30.], [':', '-', '--']):
        plt.semilogx(cycles, hettler(cycles, load), line + 'm', lw=2)
        plt.semilogx(cycles, shenton_2(cycles, load), line + 'c', lw=2)
        frequencies = [5., 10.]
        colors = ['r', 'b', 'g', 'k']
        symbols = ['x', 'o', 's', 'd']
        for f, c, sym in zip(frequencies, colors, symbols):
            settlement = [0, 0]
            for k, n in enumerate([1e3, 1e6]):
                step_name = 'cycles_' + str(int(n))
                odb_filename = (odb_directory + '/results_sleepers_high_'
                                + str(load).replace('.', '_') + 't_' + str(int(f)) + 'Hz.odb')
                up = abq.get_data_from_path(odb_filename, path_points, 'UP', 'UP2', output_position='NODAL',
                                            step_name=step_name)
                settlement[k] = -up[0]*1000
            plt.semilogx(1e6, settlement[1] - 0*settlement[0],  c + sym, lw=3, ms=12, mew=2)

    plt.xlabel('Load Cycles [-]')
    plt.ylabel('Settlement [mm]')
    fig.set_size_inches(12., 6., forward=True)
    plt.xlim(1, 2e6)
    plt.ylim(0, 35.)

    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.07, 0.12, 0.55, box.height])
    load_labels = [
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{Loads}')[0],
        plt.plot([1., 1.], [-1, -2], ':k', lw=2, label='17.5 t')[0],
        plt.plot([1., 1.], [-1, -2], '-k', lw=2, label='22.5 t')[0],
        plt.plot([1., 1.], [-1, -2], '--k', lw=2, label='30.0 t')[0]
    ]
    legend = ax.legend(handles=load_labels, loc='upper left', bbox_to_anchor=(0.03, 1.), numpoints=1)
    plt.gca().add_artist(legend)

    labels = [
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{Frequencies}')[0],
        plt.plot([1., 1.], [-1, -2], 'rx', lw=3, label='5 Hz', ms=12, mew=2)[0],
        plt.plot([1., 1.], [-1, -2], 'bo', lw=3, label='10 Hz', ms=12, mew=2)[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=3, label='white', ms=12, mew=2)[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=3, label=r'\textbf{Empirical models}', ms=12, mew=2)[0],
        plt.plot([1., 1.], [-1, -2], 'm', lw=2, label=r'Hettler (1984)', alpha=0.5)[0],
        plt.plot([1., 1.], [-1, -2], 'k', lw=2, label=r'Selig and Waters (1994)', alpha=0.5)[0],
        plt.plot([1., 1.], [-1, -2], 'c', lw=2, label=r'Shenton (1985)', alpha=0.5)[0],
        plt.plot([1., 1.], [-1, -2], 'g', lw=2, label=r'Thom and Oakley (2006))', alpha=0.5)[0],

    ]

    legend = ax.legend(handles=labels, loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
    legend.get_texts()[3].set_color("white")
    plt.gca().add_artist(legend)
    plt.savefig(figure_directory + '/empirical_models.png', dpi=600)
    plt.show()


if __name__ == '__main__':
    main()
