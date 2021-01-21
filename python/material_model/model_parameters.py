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
                 1.25454110e+01, 1.94706613e+02, -1.01551720e-01,  8.71897483e-17,
                 7.55306376e+03, 5.14233859e+00, -2.59408513e+00]

parameters[10] = [1.81696447e+00, 3.67922225e-07, 5.25068390e-02, 7.49368438e-04,
                  1.39169448e+01, 1.86223864e+02, -3.35816150e-02,  1.33050135e-14,
                  2.38379173e+03, 5.13475409e+00, -4.00005556e-01]

parameters[20] = [1.41053055e+00, 3.32379026e-06, 8.67794008e-01, 4.52482073e-04,
                  2.00671901e+00, 1.65727190e+02, 3.06103945e-02,  7.74980965e-08,
                  4.31404416e+03,  1.33261204e+00, -1.81588568e+00]

parameters[40] = [1.40956811e+00, 1.10662626e-06, 7.38373606e-02, 4.52749231e-04,
                  1.25748667e+00, 4.76597957e+02, 2.29798533e-17,  3.09954288e-09,
                  3.65065586e+03,  1.72090903e+00, -2.84897674e+00]

parameters_common = [1.69673393e+00, 2.74160533e-07, 3.77582169e-04, 8.94507739e-04,
                     1.35766134e+01, 2.42346547e+02, 8.88304920e-01, 5.69321089e-01,
                     3.82273550e-01, -9.53054047e-03, -4.97674517e-07,  3.83036719e+03,
                     7.14264717e-01, -6.85558953e-04, 1.,  1.,
                     1.,   -1.51989904e-03,  4.14039628e-04]


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
