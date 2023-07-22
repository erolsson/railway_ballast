import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.gridspec as gridspec

from multiprocesser.multiprocesser import multi_processer

from experimental_results import fouled
from material_model import MaterialModel
from model_parameters import get_parameters
from experimental_results import sun_et_al_16

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

figure_directory = os.path.expanduser('~/railway_ballast/Figures/')

def main():
    fig = plt.figure(0)
    frequencies = [5, 10, 20, 40]
    pressures = [30, 60]
    lines = ['--', '-']
    sym = ['s', 'o']
    colors = ['r', 'b', 'g', 'k']
    deviatoric_stresses = np.linspace(1, 500, 100)
    n = np.exp(np.linspace(0, np.log(5e5)))

    for p, line, sym in zip(pressures, lines, sym):
        for c, f in zip(colors, frequencies):
            static_stress = -p*np.array([1, 1, 1, 0, 0, 0])
            strain = 0*deviatoric_stresses
            for i, q in enumerate(deviatoric_stresses):
                cyclic_stress = -q*np.array([1, 0, 0, 0, 0, 0])
                parameters = get_parameters(f, common=False)
                model = MaterialModel(parameters)
                model.update(n, cyclic_stress, static_stress)
                strain[i] = -model.deviatoric_strain()[-1, 0]

            plt.plot(deviatoric_stresses, strain, line + c, lw=2)
            data_sets = sun_et_al_16.get_data(p=p, f=f)
            for experiment in data_sets:
                if experiment.cycles[-1] > 4e5:
                    plt.plot(experiment.q, experiment.deviatoric_axial_strain()[-1], c + sym, ms=12, mew=2)
    plt.ylim(0, 0.3)
    labels = [
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{p=30 kPa}')[0],
        plt.plot([1., 1.], [-1, -2], 'ks', ms=12, label=r'Exp.')[0],
        plt.plot([1., 1.], [-1, -2], '--k', lw=2, label=r'Model')[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=2, label=r'white')[0],
        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{Frequencies}')[0],
        plt.plot([1., 1.], [-1, -2], 'r', lw=2, label=r'$f=5$ Hz')[0],
        plt.plot([1., 1.], [-1, -2], 'b', lw=2, label=r'$f=10$ Hz')[0],

        plt.plot([1., 1.], [-1, -2], 'w', label=r'\textbf{p=60 kPa}')[0],
        plt.plot([1., 1.], [-1, -2], 'ko', ms=12, label=r'Exp.')[0],
        plt.plot([1., 1.], [-1, -2], '-k', lw=2, label=r'Model')[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=2, label=r'white')[0],
        plt.plot([1., 1.], [-1, -2], 'w', lw=2, label=r'white')[0],
        plt.plot([1., 1.], [-1, -2], 'g', lw=2, label=r'$f=20$ Hz')[0],
        plt.plot([1., 1.], [-1, -2], 'k', lw=2, label=r'$f=40$ Hz')[0]
    ]

    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([0.08, 0.12, 0.45, box.height])
    plt.text(0.05, 0.92, r'\bf(a)  $\boldsymbol{N=' + '500 000' + '}$', transform=ax.transAxes)

    fig.set_size_inches(16., 6.25, forward=True)
    plt.xlabel(r'Deviatoric von Mises stress, $q$ [kPa]', fontsize=24)
    plt.ylabel('Deviatoric axial strain', fontsize=24)

    legend = ax.legend(handles=labels, loc='upper left', bbox_to_anchor=(1., 1.035), numpoints=1, ncol=2)
    legend.get_texts()[3].set_color("white")
    legend.get_texts()[10].set_color("white")
    legend.get_texts()[11].set_color("white")
    plt.gca().add_artist(legend)
    plt.savefig(figure_directory + '/deviatoric_stress_strain.png')


    plt.show()


if __name__ == '__main__':
    main()
