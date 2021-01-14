from collections import OrderedDict

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

parameters = OrderedDict()


parameters[5] = [2.78174304e+00, 5.67424219e-09, 4.79444140e-02, 1.21521435e-03,
                 1.25454110e+01, 1.94706613e+02, 0.00189335, -0.00511544,  0.0013497,   0.00132903 ]

parameters[10] = [1.81696447e+00, 3.67922225e-07, 5.25068390e-02, 7.49368438e-04,
                  1.39169448e+01, 1.86223864e+02]

parameters[20] = [1.41053055e+00, 3.32379026e-06, 8.67794008e-01, 4.52482073e-04,
                  2.00671901e+00, 1.65727190e+02]

parameters[40] = [1.40956811e+00, 1.10662626e-06, 7.38373606e-02, 4.52749231e-04,
                  1.25748667e+00, 4.76597957e+02]

parameters_common = [1.69673393e+00, 2.74160533e-07, 3.77582169e-04, 8.94507739e-04,
                     1.35766134e+01, 2.42346547e+02, 8.88304920e-01, 5.69321089e-01,
                     3.82273550e-01]


def main():
    pressure = np.linspace(0, 100, 1000)
    ep = np.linspace(0, 0.3, 1000)
    for f, par in parameters.items():
        plt.figure(0)
        plt.plot(pressure, (np.sqrt(1 + par[2]*pressure + par[3]*pressure**2)))

        plt.figure(1)
        plt.plot(ep,  par[5]*(1 - np.exp(-par[4]*ep)))
        plt.figure(2)
        plt.plot(pressure, par[1]*(230/np.sqrt(1 + par[2]*pressure + par[3]*pressure**2))**par[0])

    plt.figure(2)
    plt.xlim(20, 70)
    plt.ylim(0, 0.02)
    plt.show()


if __name__ == '__main__':
    main()
