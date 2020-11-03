import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results.experimental_results import sun_et_al_16
from material_model import MaterialModel

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    par = np.array([1.78495452e+00,  5.86282720e-07,  1.53109369e+00,  1.35530067e-02,
                    1.78556464e+00,  2.04826980e+02,  1.22575778e+00,  1.64928270e+00,
                    1, -3.38200126e-02, -3.37218691e-03, -2.00237193e-06])
    for i, f in enumerate([5., 10., 20., 40.]):
        experimental_data = sun_et_al_16.get_data(f=f)
        plt.figure(i)
        for experiment in experimental_data:
            edev = (experiment.deviatoric_axial_strain() - experiment.deviatoric_axial_strain()[0])
            ea = experiment.axial_strain - experiment.axial_strain[0]
            ev = experiment.volumetric_strain - experiment.volumetric_strain[0]
            dea = np.diff(ea)/np.diff(experiment.cycles)
            dev = np.diff(ev)/np.diff(experiment.cycles)
            dedev = np.diff(edev)/np.diff(experiment.cycles)
            plt.plot(ea, ev/ea, lw=2)

    plt.show()


if __name__ == '__main__':
    main()
