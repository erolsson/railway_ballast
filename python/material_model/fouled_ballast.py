import numpy as np

from scipy.optimize import fmin

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.gridspec as gridspec

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

volumetric_parameters = {
    0: [  2.94920716e-02,  2.31123738e+01, -6.42883747e+03, -2.18051714e+00],
    10: [ 3.56088428e-02,  2.88471642e+01, -4.12606288e+03, -3.63839679e+00],
    25: [ 4.92147816e-02,  4.45747934e+01, -2.20096875e+03, -7.43197167e+00],
    50: [ 6.71514432e-02,  1.76670361e+02, -5.31103869e+02, -3.74425178e+01],
    80: [ 6.74990739e-02,  7.19437043e+02, -1.17363760e+02, -1.59589312e+02]
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
    e0 = -experiment.volumetric_strain[0]
    e_exp = -experiment.volumetric_strain - e0
    static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
    cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
    model = MaterialModel(material_parameters=par)
    model.update(experiment.cycles, cyclic_stress, static_stress)
    model_e = -model.volumetric_strain()
    diff = (model_e - e_exp)**2
    max_idx = np.argmax(diff)
    r = np.sqrt(np.sum(diff))
    return r


def determine_axial_parameters(experiment):
    return fmin(axial_residual, [2.17731700e+00, 3.49473023e-08, 4.57799168e+00, 2.54551805e+02],
                args=(experiment,), maxfun=1e6, maxiter=1e6)

def determine_volumetric_parameters(experiment, old_parameters, start_guess=None):
    if start_guess is None:
        start_guess = [1.98002052e-02,  1.58308480e+01, -3.50874250e+03, -1.22550803e+00]
    return fmin(volumetric_residual, start_guess,
                args=(old_parameters, experiment), maxfun=1e6, maxiter=1e6)


def main():
    colors = 'rbgkm'
    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))
    plt.figure(0, figsize=(14, 7))

    gs = gridspec.GridSpec(3, 2)
    parameter_job_list = []
    ax1 = plt.subplot(gs[0:2, 0:1])
    plt.xlabel('Cycles', fontsize=24)
    plt.ylabel('Deviatoric axial strain', fontsize=24)
    plt.xlim(1, 1e6)
    plt.ylim(0, 0.3)
    plt.tight_layout()

    ax2 = plt.subplot(gs[0:2, 1:2])
    plt.xlabel('Cycles', fontsize=24)
    plt.ylabel('Volumetric strain', fontsize=24)
    plt.xlim(1, 1e6)
    plt.ylim(-0.015, 0.01)
    ax2.yaxis.set_label_coords(-0.12, 0.5)
    plt.tight_layout()

    for (vci, experiment), c in zip(fouled.items(), colors):
        par = assign_axial_parameters(axial_parameters[vci])
        par = assign_volumetric_parameters(par, volumetric_parameters[vci])
        static_stress = -experiment.p*np.array([1, 1, 1, 0, 0, 0])
        cyclic_stress = -experiment.q*np.array([1, 0, 0, 0, 0, 0])
        model_1 = MaterialModel(par)
        model_1.update(cycles, cyclic_stress, static_stress)
        ea_1 = -model_1.deviatoric_strain()[:, 0]

        ax1.semilogx(experiment.cycles, experiment.deviatoric_axial_strain(), '-' + c, lw=2)

        ax1.semilogx(cycles, ea_1 + experiment.deviatoric_axial_strain()[0],
                     '--' + c, lw=2)

        ax2.semilogx(experiment.cycles, -experiment.volumetric_strain, '-' + c, lw=2)
        ax2.semilogx(cycles, -model_1.volumetric_strain() - experiment.volumetric_strain[0],
                     '--' + c, lw=2)
    lines = [plt.plot([0, -1], [-1, -1], 'w', lw=2, label=r'\textbf{VCI}')[0]]
    for vci, c in zip(fouled, colors):
        lines.append(plt.plot([0, -1], [-1, -1], c, lw=2, label=str(int(vci)) + r'\%')[0])
    plt.legend(handles=lines, ncol=2, bbox_to_anchor=(-0.8, -0.2), loc='upper left')
    plt.savefig('../../Figures/fouled_ballast.png', dpi=600)

    plt.show()


if __name__ == '__main__':
    main()