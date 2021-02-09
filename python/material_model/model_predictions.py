import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results import sun_et_al_16

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    base_parameters = np.ones(19)
    cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))
    colors = {(10, 230): 'b', (30, 230): 'r', (60, 230): 'g', (60, 370): 'k', (30, 276): 'm', (60, 460): 'y'}
    plt.figure(0, figsize=(12, 12))
    fig, axes1 = plt.subplots(2, 2, num=0)
    # axes1 = [ax1, ax2, ax3, ax4]
    # ax.set_position([0.1, 0.12, 0.55, box.height])
    for i, f in enumerate([5., 10., 20., 40.]):
        experimental_data = sun_et_al_16.get_data(f=f)
        fig_idx = i//2, i % 2
        ax = axes1[fig_idx[0]][fig_idx[1]]
        for experiment in experimental_data:
            p = experiment.p
            q = experiment.q
            ax.semilogx(experiment.cycles, experiment.deviatoric_axial_strain(),
                        '-' + colors[(p, q)], lw=2)
            plt.ylim(0, 0.3)
    plt.figure(0)
    for _, c in colors.items():
        plt.plot([-2, -1], [-1, -1], c, lw=2, label='test')
    fig.legend()
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()
