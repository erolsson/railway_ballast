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
    res = 0
    for sim in simulations:
        r, e_sim, e_exp, p, q, f = sim
        res += r
        print("\tf = " + str(f) + ",\t p = " + str(p) + " , q = " + str(q) + ": e_exp = " + str(e_exp) + "\t e_sim = "
              + str(e_sim) + "\t r = " + str(r))
    print(fitting_parameters)
    print(res)
    return res


def calc_deviatoric_residual_for_data_experiment(par, experiment):
    e0 = experiment.deviatoric_axial_strain()[0]
    e_exp = experiment.deviatoric_axial_strain() - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.deviatoric_strain()[:, 0]
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) != 0)
    r = np.sqrt(np.sum((1 - model_e[idx]/e_exp[idx])**2*np.log(experiment.cycles[idx])/e_exp.shape[0])/e_exp.shape[0])
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def calc_volumetric_residual_for_data_experiment(par, experiment):
    e0 = experiment.volumetric_strain[0]
    e_exp = experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par, frequency=experiment.f)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) > 1e-3)
    r = np.sqrt(np.sum((model_e[idx] - e_exp[idx])**2/e_exp.shape[0]))*100
    # if experiment.p == 10 and experiment.f == 10 and experiment.q == 230:
    #     r = 0
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def main():
    f = 5
    frequencies = [5, 10, 20, 40]
    # parameters_to_fit = list(range(4, 6))
    # parameters_to_fit = list(range(9, 14)) + list(range(17, 19))
    # parameters_to_fit = list(range(9))
    parameters_to_fit = list(range(9, 20))
    # parameters_to_fit = [9]
    # parameters_to_fit = range(6)

    fitting_dataset = sun_et_al_16.get_data(f=frequencies)
    par = np.array(parameters_common)
    # par = np.zeros(20)
    # par[14:17] = 0
    # par[0:9] = parameters_common[0:9]
    # par[0:6] = parameters[f][0:6]
    # par[6:9] = 1.
    # par[10] = 13
    # par[14:17] = 1.
    # par[9:14] = parameters[20][6:11]
    # par[17:19] = parameters[20][11:13]
    print(par)
    for i in range(50):
        par[parameters_to_fit] = fmin(residual, [par[parameters_to_fit]],
                                      args=(par, parameters_to_fit, fitting_dataset,
                                            calc_volumetric_residual_for_data_experiment), maxfun=1e6,
                                      maxiter=1e6)
        print(par[parameters_to_fit])


if __name__ == '__main__':
    main()
