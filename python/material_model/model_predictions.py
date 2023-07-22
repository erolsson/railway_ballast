import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.style

from experimental_results import sun_et_al_16
from material_model import MaterialModel
from model_parameters import get_parameters

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    base_parameters = np.zeros(20)
    base_parameters[6:9] = 1.
    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))
    colors = {(10, 230): 'b', (30, 230): 'r', (30, 276): 'm', (60, 230): 'g',  (60, 370): 'k', (60, 460): 'y'}
    plt.figure(0, figsize=(14, 12))
    plt.figure(1, figsize=(14, 12))
    gs = gridspec.GridSpec(14, 2)
    y_lim_vol = [(-0.03, 0.04), (-0.02, 0.05), (-0.01, 0.06), (0., 0.16)]
    for i, f in enumerate([5., 10., 20., 40.]):
        experimental_data = sun_et_al_16.get_data(f=f)
        fig_idx = (i//2, i % 2)
        par1 = get_parameters(frequency=f)
        plt.figure(0)
        ax1 = plt.subplot(gs[5*fig_idx[0] + fig_idx[0]:5*fig_idx[0] + 5 + fig_idx[0], fig_idx[1]:fig_idx[1] + 1])
        ax1.text(2, 0.27, r'\bf{$\boldsymbol  f$=' + str(int(f)) + ' Hz}')
        plt.ylim(0, 0.3)
        plt.xlabel('Cycles')
        plt.ylabel('Deviatoric axial strain')
        plt.figure(1)
        ax2 = plt.subplot(gs[5*fig_idx[0] + fig_idx[0]:5*fig_idx[0] + 5 + fig_idx[0], fig_idx[1]:fig_idx[1] + 1])
        # ax2.text(2, y_lim_vol[i][0] + 0.27/0.30*(y_lim_vol[i][1] - y_lim_vol[i][0]),
        #          r'\bf{$\boldsymbol f$=' + str(int(f)) + ' Hz}')
        ax2.text(2, -0.03 + 0.27*0.19/0.3, r'\bf{$\boldsymbol f$=' + str(int(f)) + ' Hz}')
        plt.xlabel('Cycles')
        plt.ylabel('Volumetric strain')
        # plt.ylim(*y_lim_vol[i])
        plt.ylim(-0.03, 0.16)

        for experiment in experimental_data:
            p = experiment.p
            q = experiment.q
            static_stress = -p*np.array([1, 1, 1, 0, 0, 0])
            cyclic_stress = -q*np.array([1, 0, 0, 0, 0, 0])
            model_1 = MaterialModel(par1)
            model_1.update(cycles, cyclic_stress, static_stress)

            ea_1 = -model_1.deviatoric_strain()[:, 0]

            ax1.semilogx(experiment.cycles, experiment.deviatoric_axial_strain(),
                         '-' + colors[(p, q)], lw=2)

            ax1.semilogx(cycles, ea_1 + experiment.deviatoric_axial_strain()[0],
                         '--' + colors[(p, q)], lw=2)

            ax2.semilogx(experiment.cycles, experiment.volumetric_strain,
                         '-' + colors[(p, q)], lw=2)

            ax2.semilogx(cycles[abs(ea_1) < 0.3],
                         -model_1.volumetric_strain()[abs(ea_1) < 0.3] + experiment.volumetric_strain[0],
                         '--' + colors[(p, q)], lw=2)

            ax2.yaxis.set_label_coords(-0.12, 0.5)

    for fig in [0, 1]:
        plt.figure(fig)
        for (p, q), c in colors.items():
            plt.plot([-2, -1], [-1, -1], c, lw=2, label='$p_s={p}$ kPa, $q={q}$ kPa'.format(p=p, q=q))
        plt.plot([-2, -1], [-1, -1], 'k', lw=2, label='Experiment')
        plt.plot([-2, -1], [-1, -1], '--k', lw=2, label='Model')
        plt.legend(ncol=3, bbox_to_anchor=(-1.35, -0.2), loc='upper left', columnspacing=0.7, handletextpad=0.5)
    plt.figure(0)
    plt.savefig('../../Figures/axial_strain.png', dpi=600)
    plt.figure(1)
    plt.savefig('../../Figures/volumetric_strain.png', dpi=600)
    plt.show()


if __name__ == '__main__':
    main()
