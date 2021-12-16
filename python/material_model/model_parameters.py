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
frequency_levels = [5, 10, 20, 40]

parameters[5] = [2.78174304e+00, 5.67424219e-09, 4.79444140e-02, 1.21521435e-03,
                 1.25454110e+01, 1.94706613e+02, 1.46936052e-01,  4.51955283e+01,
                 -5.12030540e+03, -7.15689757e+00, -1.27848215e+00,  8.65098143e-03,
                 -1.06993520e-04]

parameters[10] = [1.81696447e+00, 3.67922225e-07, 5.25068390e-02, 7.49368438e-04,
                  1.39169448e+01, 1.86223864e+02,  7.73461328e-02,  4.71392820e+01,
                  -9.82332826e+03, -6.38557798e+00, -5.04045480e+00,  6.22446188e-03,
                  -8.26719358e-05]

parameters[20] = [2.17731700e+00, 3.49473023e-08, 3.817471e-02, 9.20402355e-04,
                  4.57799168e+00, 2.54551805e+02, 1.98002052e-02,  1.58308480e+01,
                  -3.50874250e+03, -1.22550803e+00,  8.75142656e-01,  1.38206042e-03,
                  -1.56292371e-05]

parameters[40] = [1.40956811e+00, 1.10662626e-06, 7.38373606e-02, 4.52749231e-04,
                  1.25748667e+00, 4.76597957e+02, -8.31939736e-03,  2.75287534e+01,
                  -1.37410713e+03, -3.70763418e+00,  1.70852044e+00,  1.61342945e-03,
                  -5.94414120e-05]


def get_parameters(frequency=None, velocity=None, axle_distance=3., common=False):
    return parameters[frequency]


def get_parameters_for_frequency(f):
    par = np.zeros(19)
    par[6:9] = 1
    par[0:6] = parameters[f][0:6]
    par[9:14] = parameters[f][6:11]
    par[17:19] = parameters[f][11:13]
    return par


def main():
    p = 10
    q = 230
    ev = np.linspace(0, 0.05, 1000)
    c = 'brgk'
    for i, (f, par) in enumerate(parameters.items()):
        if f < 20:
            plt.plot(ev, par[7]*((10 + q)/2 + par[10]*q - par[8]*ev)**par[9], c[i], lw=2)
            plt.plot(ev, par[7]*((30 + q)/2 + par[10]*q - par[8]*ev)**par[9], '--' + c[i], lw=2)
            plt.plot(ev, par[7]*((60 + q)/2 + par[10]*q - par[8]*ev)**par[9], ':' + c[i], lw=2)
    plt.show()


if __name__ == '__main__':
    main()
