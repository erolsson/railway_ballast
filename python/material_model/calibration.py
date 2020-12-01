import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import fmin

from experimental_data import sun_et_al_16
from material_model import MaterialModel
from multiprocesser.multiprocesser import multi_processer


matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def deviatoric_residual(fitting_parameters, parameters, parameters_to_fit, experiments, residual_func):
    parameters[parameters_to_fit] = fitting_parameters
    job_list = [(residual_func, [parameters, experiment], {})
                for experiment in experiments]
    simulations = multi_processer(job_list, delay=0.)
    residual = 0
    for sim in simulations:
        r, e_sim, e_exp, p, q, f = sim
        residual += r
        print("\tf = " + str(f) + ",\t p = " + str(p) + " , q = " + str(q) + ": e_exp = " + str(e_exp) + "\t e_sim = "
              + str(e_sim) + "\t r = " + str(r))
    print(fitting_parameters)
    print(residual)
    return residual


def calc_deviatoric_residual_for_data_experiment(parameters, experiment):
    e0 = experiment.deviatoric_axial_strain()[0]
    e_exp = experiment.deviatoric_axial_strain() - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(parameters=parameters, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.deviatoric_strain()[:, 0]
    r = np.sum((1 - model_e[e_exp != 0]/e_exp[e_exp != 0])**2*np.log(experiment.cycles[e_exp != 0]))
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def calc_volumetric_residual_for_data_experiment(parameters, experiment):
    e0 = experiment.volumetric_strain[0]
    e_exp = experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(parameters=parameters, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    idx = np.logical_and(e_exp < 0.5, e_exp != 0)
    r = np.sum((model_e[idx] - e_exp[idx])**2*np.log(experiment.cycles[idx]))

    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def main():
    frequencies = [10.]
    parameters_to_fit = range(9, 13)
    parameters = np.array([1.92595851e+00, 1.33773133e-07, 2.36533127e-02, 6.54720828e-04,
                           1.59528891e+01, 2.03976492e+02, 1., 1.,
                           1., 2.21751785e-01, -4.55843109e-03, 1, 1])

    fitting_dataset = sun_et_al_16.get_data(f=frequencies)
    for i in range(10):
        parameters[parameters_to_fit] = fmin(deviatoric_residual, [parameters[parameters_to_fit]],
                                             args=(parameters, parameters_to_fit, fitting_dataset,
                                                   calc_volumetric_residual_for_data_experiment), maxfun=1e6,
                                             maxiter=1e6)
        print(parameters[parameters_to_fit])


if __name__ == '__main__':
    main()
