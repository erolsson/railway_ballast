from __future__ import print_function, division

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
    frequencies = [5, 10, 20, 40]

    lines = ['-', '--']
    colors = ['r', 'b', 'g', 'k']
    cycles = 100000
    pressures = np.linspace(0, 30, 100)
    strain = 0.05
    n = np.exp(np.linspace(0, np.log(cycles)))
    for j, f in enumerate(frequencies):
        parameters = [get_parameters(f, common=False), get_parameters(f, common=True)]
        for i, par in enumerate(parameters):
            job_list = []
            for p in pressures:
                job_list.append((calc_deviator, [strain, n, p, par],
                                 {}))
            q = np.array(multi_processer(job_list, delay=0., timeout=3600, cpus=12))
            plt.plot(pressures, q, lines[i] + colors[j], lw=2)
    plt.xlabel('$I_{1, static}$ [kPa]')
    plt.ylabel(r'$\sigma_{vM, cyclc}$ [kPa]')
    plt.tight_layout()
    # plt.savefig('strain_life.png')
    plt.show()


if __name__ == '__main__':
    main()
