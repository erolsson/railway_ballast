import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results.experimental_results import sun_et_al_16

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    data = sun_et_al_16.get_data(f=[40.])
    for data_set in data:
        e_dev = data_set.deviatoric_axial_strain() - data_set.deviatoric_axial_strain()[0]
        plt.figure(1)
        plt.semilogx(data_set.cycles, e_dev, lw=2)
        plt.figure(2)
        plt.plot(data_set.cycles[data_set.cycles < 1000], e_dev[data_set.cycles < 1000])
        plt.figure(3)
        dedn = np.diff(e_dev)/np.diff(data_set.cycles)
        plt.loglog(e_dev[1:], dedn)
    plt.show()


if __name__ == '__main__':
    main()

