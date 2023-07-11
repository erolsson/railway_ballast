import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib.style

from multiprocesser.multiprocesser import multi_processer

from experimental_results import fouled
from material_model import MaterialModel
from model_parameters import get_parameters

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


axial_parameters = {
    0: [0.7007383168447257, 7.88563843818382e-05, 80.97644241537431, 189.45729914239723],
    10: [0.7874659336720891, 0.00013459465181625324, 68.9902080448846, 189.4574616496355],
    25: [1.1720370643975135, 2.0150452824347696e-05, 40.33091024599454, 189.49183148736688],
    50: [4.570481724219977, 3.2445937476993676e-13, 8.069490146171887, 212.2397864327162],
    80: [1.669158786097995, 1.890051243976514e-06, 25.226468203857685, 189.4817604361287]
}

def assign_axial_parameters(new_axial_parameters):
    par = get_parameters(20.)
    par[0] = new_axial_parameters[0]
    par[1] = new_axial_parameters[1]
    par[4] = new_axial_parameters[2]
    par[5] = new_axial_parameters[3]
    return par


def assign_volumetric_parameters(old_parameters, new_volumetric_parameters):
    old_parameters[6] = new_volumetric_parameters[0]
    old_parameters[7] = new_volumetric_parameters[1]
    old_parameters[8] = new_volumetric_parameters[2]
    old_parameters[9] = new_volumetric_parameters[3]
    return old_parameters


def axial_residual(calibration_parameters, experiment):
    par = assign_axial_parameters(calibration_parameters)
    e0 = experiment.deviatoric_axial_strain()[0]
    e_exp = experiment.deviatoric_axial_strain() - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.deviatoric_strain()[:, 0]
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) != 0)
    r = np.sqrt(np.sum((1 - model_e[idx]/e_exp[idx])**2*np.log(experiment.cycles[idx])/e_exp.shape[0])/e_exp.shape[0])
    # print(r, calibration_parameters)
    return r


def volumetric_residual(calibration_parameters, old_parameters, experiment):
    par = assign_volumetric_parameters(old_parameters, calibration_parameters)
    e0 = experiment.volumetric_strain[0]
    e_exp = experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    idx = np.logical_and(e_exp < 0.3, abs(e_exp) > 1e-3)
    r = np.sqrt(np.sum((model_e[idx] - e_exp[idx])**2/e_exp.shape[0]))*100
    return r


def determine_axial_parameters(experiment):
    return fmin(axial_residual, [2.17731700e+00, 3.49473023e-08, 4.57799168e+00, 2.54551805e+02],
                args=(experiment,), maxfun=1e6, maxiter=1e6)

def determine_volumetric_parameters(experiment, old_parameters):
    return fmin(volumetric_residual, [1.98002052e-02,  1.58308480e+01, -3.50874250e+03, -1.22550803e+00],
                args=(old_parameters, experiment), maxfun=1e6, maxiter=1e6)


def main():
    colors = 'rbgkm'
    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))
    parameter_job_list = []
    for (vci, experiment), c in zip(fouled.items(), colors):
        plt.figure(0)
        plt.semilogx(experiment.cycles, experiment.deviatoric_axial_strain(), '-' + c, lw=2)

        par = assign_axial_parameters(axial_parameters[vci])
        static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
        cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
        model_1 = MaterialModel(par)
        model_1.update(cycles, cyclic_stress, static_stress)
        ea_1 = -model_1.deviatoric_strain()[:, 0]

        plt.semilogx(cycles, ea_1 + experiment.deviatoric_axial_strain()[0],
                     '--' + c, lw=2)

        plt.figure(1)
        plt.semilogx(experiment.cycles, -experiment.volumetric_strain, '-' + c, lw=2)
        plt.semilogx(cycles, -model_1.volumetric_strain() - experiment.volumetric_strain[0],
                     '--' + c, lw=2)
        parameter_job_list.append([determine_volumetric_parameters, [experiment, par], {}])
    pars = multi_processer(parameter_job_list, timeout=1e6, delay=0.)
    for par in pars:
        print(par)
    plt.show()


if __name__ == '__main__':
    main()