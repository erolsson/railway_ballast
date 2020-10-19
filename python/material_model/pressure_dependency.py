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


fig = 0
for f in [5., 10., 20., 40.]:
    experimental_data = sun_et_al_16.get_data(f=[f])
    plt.figure(fig)
    for experiment in experimental_data:
        dedn = np.diff(experiment.strain)/np.diff(experiment.cycles)

        plt.semilogx(experiment.cycles[1:], dedn, lw=2)
    fig += 1

plt.show()
