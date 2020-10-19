from __future__ import print_function, division

import sys
import warnings
import find_modules  # noqa

from common import numpy as np
from common import scipy
from multiprocesser.multiprocesser import multi_processer

from scipy.integrate import solve_ivp
from scipy.optimize import fmin


freq_idx = {10: 7, 20: 8, 40: 9}


def permanent_strain(cycles, p, q, freq, parameters, e0=0):
    scalar = False
    if not isinstance(cycles, np.ndarray):
        cycles = np.array([cycles])
        scalar = True
    strain = 0*cycles
    # parameters = abs(parameters)

    def dedn(n, ep):
        gf = parameters[0]
        A = parameters[1]
        A1 = parameters[2]
        A2 = abs(parameters[3])
        arg = 1. + A1*p + A2*p**2
        if freq == 5.:
            freq_factor = 1.
        else:
            freq_factor = parameters[freq_idx[freq]]
        if arg < 1e-6:
            arg = 1e-6
        f = (freq_factor*q/np.sqrt(arg) - hf(ep))
        # f = (q/p - hf(ep))
        e = 0*f
        e[f > 0] = A*f[f > 0]**gf
        return e

    def hf(ep):
        # nf = parameters[4]
        # H0 = parameters[5]
        # H1 = parameters[6]
        # return H0 + H1*(1-np.exp(-nf*(ep - e0)))
        f_vals = [8.13220695e-02,
  2.68868736e+01,  4.91914262e+01, -2.45166693e+01, -2.33720631e+01,
 -4.08586276e+01, -1.00691061e-02, -4.69099710e-01,  6.33116729e+01,
 -4.92231884e-01, -5.99374859e+01,  4.48312718e+00, -4.41065771e-02,
 -5.29131981e+00, -1.17381751e+02,  4.53602062e+00, -9.32487550e+00,
  9.18212675e-02,  8.05041126e+00,  1.19068324e+00]

        strain_values = np.linspace(0, 0.3, 20)
        f_vals = np.abs(f_vals)
        return np.interp(ep, strain_values, np.cumsum(f_vals))

    for i in range(cycles.shape[0]):
        if i == 0:
            n0 = 1
            e1 = e0
        else:
            e1 = strain[i-1]
            n0 = strain[i-1]

        solution = solve_ivp(dedn, [n0, cycles[i]], [e1])
        strain[i] = solution.y[0][-1]
        if strain[i] > 10:
            strain[i] = 10
    if scalar:
        strain = strain[0]
    return strain


def calc_residual_for_dataset(points, data, par):
    e0 = np.interp(np.log(1), np.log(data.cycles), data.strain)
    model_e = permanent_strain(points, data.p, data.q, data.f, par, e0=0)
    experiment_e = np.interp(np.log(points), np.log(data.cycles), data.strain)
    return np.sum((model_e - (experiment_e - e0))**2)


def residual(parameters, datasets):
    job_list = []
    for dataset in datasets:
        eval_points = np.exp(np.linspace(np.log(1), np.log(np.max(dataset.cycles)), 100))
        job_list.append((calc_residual_for_dataset, [eval_points, dataset, parameters], {}))
    residuals = multi_processer(job_list, delay=0, cpus=1)
    r = sum(residuals)/len(residuals)
    print(parameters)
    print(residuals)
    print(r)
    return r


def main():
    import matplotlib.pyplot as plt
    import matplotlib.style

    matplotlib.style.use('classic')
    plt.rc('text', usetex=True)
    plt.rc('font', serif='Computer Modern Roman')
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                      'monospace': ['Computer Modern Typewriter']})
    from experimental_results.experimental_results import sun_et_al_16
    fig = plt.figure(0)
    data = sun_et_al_16.get_data(f=[5., 10., 20., 40])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        par = fmin(residual, [4.35165364e-01,  4.88334854e-05,  1.01199265e+03,  2.00504418e-01,
  5.95614711e+00, -2.60338013e-01,  5.50206239e+00,  1.45525929e+00,
  1.53855819e+00,  1.82715745e+00],
                   args=(data,), maxfun=1e6, maxiter=1e6)
    # par = np.array([3.26067740e+00  1.64682732e-06  6.24063749e+01 -7.27014619e-13
    #   8.67609593e+00  1.35453333e+01])

    # 2.01266553e+00  1.06441488e-06  7.78057288e+00 -2.40037526e-13
    #   8.17380463e+00  3.32964592e+01

    # 2.18980083e+00  2.61294997e-06  1.71552211e+01 -1.04617517e-12
    #   1.27409840e+00  6.92076296e+01

    # 1.23019992e+00  1.29236941e-04  2.30681914e+04 -6.42318544e-10
    # 7.33674237e-04  5.28999887e+03  1.23997407e+00  1.80326425e+00
    # 2.77779158e+00
    #
    print(par)
    N = np.exp(np.linspace(0, np.log(1e6)))
    color = ['r', 'b', 'g', 'k']
    for experiment, c in zip(data, color):
        print(experiment.p, experiment.q)
        model_strain = permanent_strain(experiment.cycles, experiment.p, experiment.q, par)
        if c == 'k':
            plt.semilogx(experiment.cycles, experiment.strain, '-' + c, lw=2, label='Experiment (Sun et. al 2016)')
            plt.semilogx(experiment.cycles, model_strain + experiment.strain[0], '--' + c, lw=2, label='Model')
        else:
            plt.semilogx(experiment.cycles, experiment.strain, '-' + c, lw=2)
            plt.semilogx(experiment.cycles, model_strain + experiment.strain[0], '--' + c, lw=2)
    # plt.ylim(0, 0.15)
    plt.xlim(0, 5e5)
    fig.set_size_inches(10, 8., forward=True)
    plt.legend(loc='best')
    plt.text(3e3, 0.02, '$p=60$ kPa, $q=230$ kPa')
    plt.text(3e3, 0.085, '$p=60$ kPa, $q=370$ kPa', color='g')
    plt.text(3e3, 0.045, '$p=30$ kPa, $q=230$ kPa', color='r')
    plt.text(1200, 0.11, '$p=10$ kPa, $q=230$ kPa', color='b')
    plt.xlabel('Cycles $N$')
    plt.ylabel(r'Permanent strain $\varepsilon_p$')
    plt.tight_layout()
    plt.savefig('../../Figures/permanent_strain_10_60.png')
    plt.show()


if __name__ == '__main__':
    main()
