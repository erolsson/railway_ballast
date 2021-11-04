import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from comparison_of_models import get_path_points_for_fem_simulation
from abaqus_python_interface import ABQInterface

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 24})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')

abq = ABQInterface("abq2018")


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
    cycles = np.array([float(10**i) for i in range(7)])
    plt.semilogx(cycles, thom_1(cycles, 0), ':k', lw=2)
    plt.semilogx(cycles, selig_2(cycles, 0), ':k', lw=2)
    path_points = get_path_points_for_fem_simulation('sleepers_high')

    for load, line in zip([17.5], [':', '-', '--']):
        plt.semilogx(cycles, hettler(cycles, load), line + 'k', lw=2)
        plt.semilogx(cycles, shenton_2(cycles, load), line + 'k', lw=2)

    for rail_fixture, c in zip(['slab', 'sleepers'], ['b', 'r']):
        odb_filename = (odb_directory + '/results_' + rail_fixture + '_high_'
                        + str(17.5).replace('.', '_') + 't_' + str(int(5)) + 'Hz.odb')
        settlement = 0*cycles
        for k, n in enumerate(cycles[1:], 1):
            step_name = 'cycles_' + str(int(n))
            up = abq.get_data_from_path(odb_filename, path_points, 'UP', 'UP2', output_position='NODAL',
                                        step_name=step_name)
            settlement[k] = -up[0]*1000
        plt.semilogx(cycles[settlement != np.nan], settlement[settlement != np.nan], '-' + c,
                     lw=2, ms=12, mew=2, label=rail_fixture.capitalize())

    plt.semilogx([1, 2], [-1, -1], ':k', lw=2, label="Empirical models")
    plt.xlabel('Load Cycles [-]')
    plt.ylabel('Settlement [mm]')
    fig.set_size_inches(12., 8., forward=True)
    plt.xlim(1, 2e6)
    plt.ylim(0, 35.)
    plt.legend(loc='best')
    plt.savefig(figure_directory + '/empirical_models.png')
    plt.show()


if __name__ == '__main__':
    main()
