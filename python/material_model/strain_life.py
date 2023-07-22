from __future__ import print_function, division

import os

import matplotlib.pyplot as plt
import matplotlib.style

import numpy as np
from scipy.optimize import fmin

from material_model import MaterialModel
from model_parameters import get_parameters
from multiprocesser.multiprocesser import multi_processer

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

figure_directory = os.path.expanduser('~/railway_ballast/Figures/')


def residual(deviator, strain, cycles, pressure, parameters):
    static_stress = -pressure*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -deviator*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(parameters)
    model.update(cycles, cyclic_stress, static_stress)
    e = -model.deviatoric_strain()[-1, 0]
    return (e-strain)**2


def calc_deviator(strain, n, p, par):
    return fmin(residual, 10., args=(strain, n, p, par))


def main():
    fig = plt.figure(0)
    frequencies = [5, 10, 20, 40]

    lines = [':', '--', '-']
    colors = ['r', 'b', 'g', 'k']
    cycles = 5e5
    pressures = np.linspace(0, 60)
    strain_levels = [0.01, 0.05, 0.1]
    n = np.exp(np.linspace(0, np.log(cycles)))

    for j, f in enumerate(frequencies):
        par = get_parameters(f, common=False)
        for i, strain in enumerate(strain_levels):
            job_list = []
            for p in pressures:
                job_list.append((calc_deviator, [strain, n, p, par],
                                 {}))
            q = np.array(multi_processer(job_list, delay=0., timeout=3600, cpus=12))
            plt.plot(pressures, q, lines[i] + colors[j], lw=2, label='$f=' + str(f) + '$ Hz')

    labels = [
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{Strain levels}')[0],
        plt.plot([1., 1.], [-1, -2], ':k', lw=2, label=r'$\varepsilon_{eff}=0.01$')[0],
        plt.plot([1., 1.], [-1, -2], '--k', lw=2, label=r'$\varepsilon_{eff}=0.05$')[0],
        plt.plot([1., 1.], [-1, -2], '-k', lw=2, label=r'$\varepsilon_{eff}=0.10$')[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=2, label=r'white')[0],
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{Frequencies}')[0],
        plt.plot([1., 1.], [-1, -2], 'r', lw=2, label=r'$f=5$ Hz')[0],
        plt.plot([1., 1.], [-1, -2], 'b', lw=2, label=r'$f=10$ Hz')[0],
        plt.plot([1., 1.], [-1, -2], 'g', lw=2, label=r'$f=20$ Hz')[0],
        plt.plot([1., 1.], [-1, -2], 'k', lw=2, label=r'$f=40$ Hz')[0]
    ]
    plt.xlim(0, 60)
    plt.ylim(0)
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.1, 0.12, 0.55, box.height])

    plt.xlabel('Static pressure, $p_s$ [kPa]', fontsize=24)
    plt.ylabel(r'Deviatoric von Mises stress, $q$ [kPa]', fontsize=24)
    plt.text(0.05, 0.92, r'\bf(b)  $\boldsymbol{N=' + '500 000' + '}$', transform=ax.transAxes)
    # plt.text(0.05, 0.85, r'$\boldsymbol{N=' + str(int(cycles)) + '}$', transform=ax.transAxes)
    fig.set_size_inches(12., 6.25, forward=True)

    legend = ax.legend(handles=labels, loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1)
    legend.get_texts()[4].set_color("white")
    plt.gca().add_artist(legend)
    plt.savefig(figure_directory + '/strain_yield_surface_1e6.png')
    plt.show()


if __name__ == '__main__':
    main()
