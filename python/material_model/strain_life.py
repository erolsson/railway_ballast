from __future__ import print_function, division

import matplotlib.pyplot as plt
import matplotlib.style

import numpy as np
from scipy.optimize import fmin

from functional_shapes import permanent_strain
from multiprocesser.multiprocesser import multi_processer

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def residual(deviator, strain, cycles, pressure, parameters):
    return (strain - permanent_strain(cycles, pressure, deviator, parameters)[-1])**2


def calc_deviator(strain, n, p, par):
    return fmin(residual, 10., args=(strain, n, p, par))


def main():
    par = np.array([3.78059234e+00, 2.28236862e-06, 8.14579983e+01, -9.66505289e-13, 6.80152840e+00,
                    1.28249309e+01])

    cycles = 100000

    pressures = np.linspace(0, 60, 100)
    strain_levels = [0.01, 0.05, 0.1]
    n = np.exp(np.linspace(0, np.log(cycles)))
    for strain in strain_levels:
        job_list = []
        for i, p in enumerate(pressures):
            job_list.append((calc_deviator, [strain, n, p, par],
                             {}))
        q = np.array(multi_processer(job_list, delay=0., timeout=3600, cpus=12))
        plt.plot(-3*pressures, q, lw=2)
    plt.xlabel('$I_{1, static}$ [kPa]')
    plt.ylabel(r'$\sigma_{vM, cyclc}$ [kPa]')
    plt.tight_layout()
    plt.savefig('strain_life.png')
    plt.show()


if __name__ == '__main__':
    main()
