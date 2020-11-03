import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import fmin

from experimental_results.experimental_results import sun_et_al_16
from material_model import MaterialModel
from multiprocesser.multiprocesser import multi_processer


matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def deviatoric_residual(fitting_parameters, parameters, parameters_to_fit, experiments):
    parameters[parameters_to_fit] = fitting_parameters
    job_list = [(calc_deviatoric_residual_for_data_experiment, [parameters, experiment], {})
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
    r = np.sum((1 - model_e[e_exp > 0]/e_exp[e_exp > 0])**2*np.log(experiment.cycles[e_exp > 0]))
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def main():
    frequencies = [5, 10., 20., 40]
    parameters_to_fit = range(9)
    parameters = np.array([1.39289920e+00, 5.62356963e-07, 5.69024102e-02, 1.04363423e-03,
                           3.10392063e+00, 5.30904313e+02, 1.12082052e+00, 1.72836168e+00,
                           2.37893774e+00, -4.75186794e-02, -3.02491126e-03,  1.77430478e-06,
                           1.38965363e-06])

    fitting_dataset = sun_et_al_16.get_data(f=frequencies)
    for i in range(10):
        parameters[parameters_to_fit] = fmin(deviatoric_residual, [parameters[parameters_to_fit]],
                                             args=(parameters, parameters_to_fit, fitting_dataset), maxfun=1e6,
                                             maxiter=1e6)
        print(parameters[parameters_to_fit])


if __name__ == '__main__':
    main()
