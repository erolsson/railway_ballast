import numpy as np

from scipy.integrate import solve_ivp
from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results.experimental_results import sun_et_al_16

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def permanent_strain(cycles, p, q, parameters):
    strain = 0*cycles
    parameters = abs(parameters)

    def dedn(n, ep):
        gf = parameters[0]
        a0 = parameters[1]
        a1 = parameters[2]
        a2 = parameters[3]
        f = (q/np.sqrt(1 + a1*p + a2*p**2) - hf(ep))
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
    data = sun_et_al_16.get_data(f=5., p=[10, 30, 60])
    # par = fmin(residual, [5.67, 6.344e-7, 100, 1, 7.5, 13], args=(data, ), maxfun=1e6, maxiter=1e6)
    par = np.array([5.19502802e+00, 7.64481869e-07, 4.70091315e+01, 1.34948445e-13,
                    1.14929617e+01, 1.12439205e+01])
    N = np.exp(np.linspace(0, np.log(1e6)))
    for experiment in data:
        plt.semilogx(experiment.cycles, experiment.strain, '-', lw=2)
        model_strain = permanent_strain(experiment.cycles, experiment.p, experiment.q, par)
        plt.semilogx(experiment.cycles, model_strain, '--', lw=2)

    model_strain = permanent_strain(experiment.cycles, 10, 100, par)
    plt.semilogx(experiment.cycles, model_strain, '--', lw=2)
    plt.show()
