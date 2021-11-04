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
plt.rcParams.update({'font.size': 36})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{xcolor}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')


def main():
    rail_fixture = 'sleepers'
    geometry = 'low'
    load = 22.5
    path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
    cycles = np.array([float(10**i) for i in range(7)])
    for frequency, c, sym in zip([5., 20.], ['r', 'g'], ['x', 's']):
        odb_filename = (odb_directory + '/results_' + rail_fixture + '_' + geometry + '_'
                        + str(load).replace('.', '_') + 't_' + str(int(frequency)) + 'Hz.odb')
        settlement = 0*cycles
        for k, n in enumerate(cycles[1:], 1):
            step_name = 'cycles_' + str(int(n))
            up = abq.get_data_from_path(odb_filename, path_points, 'UP', 'UP2', output_position='NODAL',
                                        step_name=step_name)
            settlement[k] = -up[0]*1000
        plt.semilogx(cycles, settlement, '-' + c + sym, lw=3, ms=12, mew=2)
    plt.xlabel('Belastningscykler')
    plt.ylabel(r'S{\"a}ttningar [mm]')
    plt.xlim(1, 1e6)
    plt.ylim(0, 50.)
    plt.tight_layout()
    plt.savefig(figure_directory + 'one_graph_settlement.png')
    plt.show()


if __name__ == '__main__':
    main()
