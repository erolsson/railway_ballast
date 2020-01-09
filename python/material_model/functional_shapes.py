import sys
import find_modules  # noqa

from common import numpy as np
from common import scipy

from scipy.integrate import solve_ivp
from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib.style

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def permanent_strain(cycles, p, q, parameters):
    scalar = False
    if not isinstance(cycles, np.ndarray):
        cycles = np.array([cycles])
        scalar = True
    strain = 0*cycles
    print strain.shape
    parameters = abs(parameters)

    def dedn(n, ep):
        gf = parameters[0]
        a0 = parameters[1]
        a1 = parameters[2]
        a2 = parameters[3]
        arg = 1 + a1*p + a2*p**2
        arg[arg < 1e-6] = 1e-6
        f = (q/np.sqrt(arg) - hf(ep))
        # f = (q/p - hf(ep))
        e = 0*f
        e[f > 0] = a0*f**gf
        return e

    def hf(ep):
        nf = parameters[4]
        a = parameters[5]
        return (1-np.exp(-nf*ep))*a

    for i in range(cycles.shape[0]):
        if i == 0:
            e0 = 0
            n0 = 0
        else:
            e0 = strain[i-1]
            n0 = strain[i-1]
        if e0 < 1:
            solution = solve_ivp(dedn, [n0, cycles[i]], [e0])
            strain[i] = solution.y[0][-1]
        else:
            strain[i] = 1.
    if scalar:
        strain = strain[0]
    return strain


def residual(parameters, *args):
    datasets = args[0]
    r = 0
    for dataset in datasets:
        model_e = permanent_strain(dataset.cycles, dataset.p, dataset.q, parameters)
        r += np.sum((dataset.strain - model_e)**2)
    print parameters, r
    return r


if __name__ == '__main__':
    from experimental_results.experimental_results import sun_et_al_16
    fig = plt.figure(0)
    data = sun_et_al_16.get_data(f=5., p=[10, 30, 60])
    # par = fmin(residual, [5.67, 6.344e-7, 100, 1, 7.5, 13], args=(data, ), maxfun=1e6, maxiter=1e6)
    # par = np.array([5.19502802e+00, 7.64481869e-07, 4.70091315e+01, 1.34948445e-13,
    #               1.14929617e+01, 1.12439205e+01])
    par = np.array([3.78059234e+00, 2.28236862e-06,  8.14579983e+01, -9.66505289e-13,   6.80152840e+00,
                    1.28249309e+01])
    N = np.exp(np.linspace(0, np.log(1e6)))
    color = ['r', 'b', 'g', 'k']
    for experiment, c in zip(data, color):
        print experiment.p, experiment.q
        model_strain = permanent_strain(experiment.cycles, experiment.p, experiment.q, par)
        if c == 'k':
            plt.semilogx(experiment.cycles, experiment.strain, '-' + c, lw=2, label='Experiment (Sun et. al 2016)')
            plt.semilogx(experiment.cycles, model_strain, '--' + c, lw=2, label='Model')
        else:
            plt.semilogx(experiment.cycles, experiment.strain, '-' + c, lw=2)
            plt.semilogx(experiment.cycles, model_strain, '--' + c, lw=2)
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
    plt.savefig('../Figures/permanent_strain_10_60.png')
    plt.show()
