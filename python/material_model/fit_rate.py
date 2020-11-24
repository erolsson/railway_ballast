from __future__ import print_function, division

import numpy as np
from scipy.optimize import fmin

from experimental_results import sun_et_al_16

import matplotlib.pyplot as plt
import matplotlib.style

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def rate_model(parameters, e_eff, p, q, frequency):
    gf = abs(parameters[0])
    A = abs(parameters[1])
    A1 = abs(parameters[2])
    A2 = parameters[3]

    nf = abs(parameters[4])
    H1 = abs(parameters[5])

    freq_idx = {10.: 6, 20.: 7, 40.: 8}
    if frequency == 5.:
        fd = 1.
    else:
        fd = parameters[freq_idx[frequency]]
    arg = 1. + A1*p + A2*p**2
    hf = H1*(1-np.exp(-nf*e_eff))
    rate = 0*hf
    rate[(fd*q/np.sqrt(arg) > hf)] = A*(fd*q/np.sqrt(arg) - hf[(fd*q/np.sqrt(arg) > hf)])**gf
    return rate


def rate_residual(parameters, experiments):
    residual = 0
    for experiment in experiments:
        e_eff = np.linspace(experiment.deviatoric_axial_strain()[0],
                            experiment.deviatoric_axial_strain()[-1], 100)
        dedn_exp = np.interp(e_eff, experiment.deviatoric_axial_strain()[1:],
                             np.diff(experiment.deviatoric_axial_strain())/np.diff(experiment.cycles))
        dedn_model = rate_model(parameters, e_eff, experiment.p,
                                experiment.q, experiment.f)
        dedn_exp[dedn_exp < 1e-12] = 1e-12
        dedn_model[dedn_model < 1e-12] = 1e-12
        residual += np.sum((np.log(dedn_model) - np.log(dedn_exp))**2)
    print(parameters)
    print(residual)
    return residual


def main():
    frequencies = [5, 10, 20, 40.]
    experiments = sun_et_al_16.get_data(f=frequencies)
    par = [3.46307406e+00,  1.61051398e-09,  2.59113736e-01,  5.78186249e-04,
           1.50087199e+01,  1.32159618e+02, 1.2,  1.5,
           1.7]

    # 5 Hz
    # 3.46307406e+00  1.61051398e-09  2.59113736e-01  5.78186249e-04
    #   1.50087199e+01  1.32159618e+02 -1.13062306e+02 -8.46062813e+00
    #   4.78849016e+02

    # 10 Hz
    # 2.63659509e+00  1.46740630e-08  9.10000445e-02  1.98652118e-04
    # 1.94577089e+01  1.68585581e+02 -2.09840120e+03 -6.93044061e+02
    # -9.68568969e+03

    # 20 Hz
    # 1.53739943e+00  1.03125336e-06  9.05338852e-02  2.10121557e-04
    # 7.17788010e+00  2.14008410e+02 -1.42881101e+05  3.16920912e+05
    #  4.04118476e+05

    for i in range(100):
        par = fmin(rate_residual, par, args=(experiments,),
                   maxfun=1e6, maxiter=1e6)
        print(par)
    for i, freq in enumerate(frequencies):
        data = sun_et_al_16.get_data(f=[freq])
        for experiment in data:
            dedn_exp = np.diff(experiment.deviatoric_axial_strain())/np.diff(experiment.cycles)
            plt.figure(i)
            plt.semilogy(experiment.deviatoric_axial_strain()[1:], dedn_exp)
            dedn_model = rate_model(par, experiment.deviatoric_axial_strain()[1:], experiment.p,
                                    experiment.q, experiment.f)
            plt.semilogy(experiment.deviatoric_axial_strain()[1:], dedn_model, '--')
    plt.show()


if __name__ == '__main__':
    main()
