import numpy as np

from scipy.optimize import fmin
from scipy.integrate import solve_ivp

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results import sun_et_al_16
from multiprocesser.multiprocesser import multi_processer

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


freq_idx = {10: 4, 20: 5, 40: 6}
num_points = 20
strain_values = np.linspace(0, 0.3, num_points)


def stab_func(par, e):
    par = np.abs(par)
    return np.interp(e, strain_values, np.cumsum(par))


def rate(par, e, p, q, f):

    if f == 5:
        freq_f = 1.
    else:
        freq_f = par[freq_idx[f]]
    g = freq_f*q/np.sqrt(1+par[2]*p + par[3]*p**2)
    arg = g-stab_func(par[7:], e)
    dedn = 0*e
    dedn[arg > 0] = par[1]*arg[arg > 0]**par[0]
    return dedn


def residual_func_(parameters, experiment):
    cycles = np.exp(np.linspace(0, np.log(experiment.cycles[-1]), 100))
    e0 = experiment.strain[0]
    exp_strain = np.interp(np.log(cycles), np.log(experiment.cycles), experiment.strain) - e0
    model_strain = permanent_strain(cycles, experiment.p, experiment.q, experiment.f, parameters)
    return np.sum((exp_strain - model_strain)**2), exp_strain[-1] + e0, model_strain[-1] + e0, experiment, cycles[-1]


def residual(par, experiments):
    job_list = []

    for experiment in experiments:
        job_list.append((residual_func_, [par, experiment], {}))
    residuals = multi_processer(job_list, delay=0.)
    r = sum([res[0] for res in residuals])
    for result in residuals:
        print("f =", result[3].f, " p =", result[3].p, " q =", result[3].q,
              "result at N =", result[4], " exp_e =", result[1], " model_e =", result[2])
    print(par)
    print(r)
    return r


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
        e = 0*f
        e[f > 0] = A*f[f > 0]**gf
        return e

    def hf(ep):
        par = np.abs(parameters)
        return np.interp(ep, strain_values, np.cumsum(par[7:]))

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


def main():
    par = [2.41831643e+00,  5.84725898e-10, -2.27076233e-03,  5.35696009e-04,
           1.13411551e+00,  1.36038327e+00,  1.87775732e+00] + [1]*num_points
    """
    par = [ 1.82921119e+00,  4.52536694e-08,  5.64026331e-02,  1.14187224e-02,
            2.07777889e+01, -6.99129677e-01,  8.97144114e-02, -1.25119288e+00,
            -3.81830553e-01, -2.18290308e+01, -2.86599225e+01,  1.23924082e+00,
            -7.82547868e+01, -2.63875799e-01, -2.30900154e+00, -5.39102083e+00,
            5.00107267e-01, -9.13250785e+00, -1.52523519e+00, -2.34323530e+00,
            2.85873951e-01,  4.40963476e+00,  4.45866027e+00, -5.15590396e+00,
            1.12170332e+01,  1.89540703e+01, -7.54540218e-02,
    ]
    """
    """
    1.86010034e+00  5.02588837e-09  5.06535457e-02  7.26990151e-04
    1.09486676e+00  1.41337306e+00  1.81160346e+00 -5.01746856e-07
    -2.07808689e-06 -9.54842974e+01  2.14716742e+00  4.14419929e+01
    -1.05189488e+01 -3.21392422e+01 -3.15880939e-07  9.18979342e-09
    -4.82989843e+01 -1.48586760e+00  2.69257251e+01 -4.45586273e-05
    -2.51294323e+01  7.36613899e+01 -2.27231699e-07 -5.04822354e-02
    -1.79516230e+00 -1.18545390e-09 -1.04609364e+02
    """
    print(', '.join([str(p) for p in par]))
    experiments = sun_et_al_16.get_data(f=[5., 10., 20., 40])

    for i in range(40):
        par = fmin(residual, par, args=(experiments,), maxfun=1e3, maxiter=1e3)

    print(', '.join([str(p) for p in par]))
    plt.figure(0)
    plt.plot(strain_values, stab_func(par[7:], strain_values))
    plt.show()


if __name__ == '__main__':
    main()
