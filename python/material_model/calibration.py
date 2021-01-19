import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import fmin

from experimental_results import sun_et_al_16
from material_model import MaterialModel
from model_parameters import parameters, parameters_common
from multiprocesser.multiprocesser import multi_processer


matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def residual(fitting_parameters, par, parameters_to_fit, experiments, residual_func):
    par[parameters_to_fit] = fitting_parameters
    job_list = [(residual_func, [par, experiment], {})
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


def calc_deviatoric_residual_for_data_experiment(par, experiment):
    e0 = experiment.deviatoric_axial_strain()[0]
    e_exp = experiment.deviatoric_axial_strain() - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.deviatoric_strain()[:, 0]
    idx = np.logical_and(e_exp < 0.3, e_exp != 0)
    r = np.sqrt(np.sum((model_e[idx] - e_exp[idx])**2)/e_exp.shape[0])*100
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def calc_volumetric_residual_for_data_experiment(par, experiment):
    e0 = experiment.volumetric_strain[0]
    e_exp = experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) > 1e-2)
    r = np.sqrt(np.sum((model_e[idx] - e_exp[idx])**2)/e_exp.shape[0])*100

    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def main():
    frequencies = [5, 10, 20, 40]
    # parameters_to_fit = list(range(4, 6))
    parameters_to_fit = list(range(9, 18))
    # parameters_to_fit = [9]
    # parameters_to_fit = range(6)
    par = np.array([1.81696447e+00, 3.67922225e-07,   5.25068390e-02, 7.49368438e-04,
                    1.39169448e+01, 1.86223864e+02,   1,              1.,
                    1,              0.08608095,  -2.48672886e-10,  2.81664104e+06,  1.09291862e+00, -9.55252021e+02,
                    1.19911441e+00,  2.12567827e+00,  3.49779516e+00, 0.])

    fitting_dataset = sun_et_al_16.get_data(f=frequencies)
    par = np.array(parameters_common)
    # par[9] = parameters[5][6]
    for i in range(6):

        par[parameters_to_fit] = fmin(residual, [par[parameters_to_fit]],
                                      args=(par, parameters_to_fit, fitting_dataset,
                                            calc_volumetric_residual_for_data_experiment), maxfun=1e6,
                                      maxiter=1e6)
        print(par[parameters_to_fit])


if __name__ == '__main__':
    main()
