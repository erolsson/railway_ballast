import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import fmin

from experimental_results import sun_et_al_16
from material_model import MaterialModel
from model_parameters import get_parameters
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
    static_stress = -30*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -320*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(np.array([1, 1e5]), cyclic_stress, static_stress)
    inst_e = -model.deviatoric_strain()[-1, 0]
    for sim in simulations:
        r, e_sim, e_exp, p, q, f = sim

        res += r
        print("\tf = " + str(f) + ",\t p = " + str(p) + " , q = " + str(q) + ": e_exp = " + str(e_exp) + "\t e_sim = "
              + str(e_sim) + "\t r = " + str(r))
    if inst_e < 0.8:
        print("Strain at problematic load", inst_e)
        res += (0.8 - inst_e)**2*100
    print(fitting_parameters)
    print(res)

    return res


def calc_deviatoric_residual_for_data_experiment(par, experiment):
    e0 = experiment.deviatoric_axial_strain()[0]
    e_exp = experiment.deviatoric_axial_strain() - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.deviatoric_strain()[:, 0]
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) != 0)
    e_corr = 0 + model_e
    e_corr[e_corr > 0.3] = 0.3
    r = np.sqrt(np.sum((e_corr[idx] - e_exp[idx])**2))
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def calc_volumetric_residual_for_data_experiment(par, experiment):
    e0 = experiment.volumetric_strain[0]
    e_exp = experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) > 1e-3)
    r = np.sqrt(np.sum((model_e[idx] - e_exp[idx])**2/e_exp.shape[0]))*100
    # if experiment.p == 10 and experiment.f == 10 and experiment.q == 230:
    #     r = 0
    return r, round(model_e[-1] + e0, 4), round(e_exp[-1] + e0, 4), experiment.p, experiment.q, experiment.f


def main():
    f = 40
    frequencies = [f]
    parameters_to_fit = [0, 1, 2, 3, 4, 5]
    fitting_dataset = sun_et_al_16.get_data(f=frequencies)
    par = np.array(get_parameters(40, common=False))

    print(par)
    for i in range(50):
        old_par = par[parameters_to_fit]
        par[parameters_to_fit] = fmin(residual, [old_par],
                                      args=(par, parameters_to_fit, fitting_dataset,
                                            calc_deviatoric_residual_for_data_experiment), maxfun=1e6,
                                      maxiter=1e6)
        if np.sum((1 - par[parameters_to_fit]/old_par)**2) < 1e-3:
            break
        print(par[parameters_to_fit])


if __name__ == '__main__':
    main()
