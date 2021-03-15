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

# parameters[20] = [1.41053055e+00, 3.32379026e-06, 8.67794008e-01, 4.52482073e-04,
#                   2.00671901e+00, 1.65727190e+02, 1.98002052e-02,  1.58308480e+01,
#                   -3.50874250e+03, -1.22550803e+00,  8.75142656e-01,  1.38206042e-03,
#                   -1.56292371e-05]


parameters[20] = [2.17731700e+00, 3.49473023e-08, 3.817471e-02, 9.20402355e-04,
                  4.57799168e+00, 2.54551805e+02, 1.98002052e-02,  1.58308480e+01,
                  -3.50874250e+03, -1.22550803e+00,  8.75142656e-01,  1.38206042e-03,
                  -1.56292371e-05]

parameters[40] = [1.40956811e+00, 1.10662626e-06, 7.38373606e-02, 4.52749231e-04,
                  1.25748667e+00, 4.76597957e+02, -8.31939736e-03,  2.75287534e+01,
                  -1.37410713e+03, -3.70763418e+00,  1.70852044e+00,  1.61342945e-03,
                  -5.94414120e-05]

parameters_common = [2.20517242e+00, 3.72924098e-08, 3.77027776e-02, 7.16206717e-04,
                    1.34590082e+01, 2.12620654e+02, 8.65661146e-01, 5.17534627e-01,
                    3.43653164e-01,  1.35287088e-01,  2.11853188e+01, -7.13235258e+03, -2.00711284e+00,
                    -2.81461850e+00,  3.82300562e-02, 1.03239331e-01,  1.32871815e-01,
                    4.19011276e-03, -4.99458815e-05]


def get_parameters(frequency=None, velocity=None, axle_distance=3., common=False):
    if frequency is not None and not common:
        try:
            return parameters[frequency]
        except KeyError:
            raise KeyError("Parameters for frequency {f} not found".format(f=frequency))
    else:
        if velocity is not None:
            frequency = velocity/3.6/axle_distance
            if frequency < frequency_levels[0] or frequency > frequency_levels[-1]:
                min_vel = frequency_levels[0]*axle_distance*3.6
                max_vel = frequency_levels[-1]*axle_distance*3.6
                if velocity < min_vel:
                    raise ValueError("The minimum allowed velocity is {v}".format(v=min_vel))
                if velocity < max_vel:
                    raise ValueError("The maximum allowed velocity is {v}".format(v=max_vel))
        par = np.zeros(13)
        par[0:6] = parameters_common[0:6]
        par[6:11] = parameters_common[9:14]
        par[11:13] = parameters_common[17:19]
        par[4] *= np.interp(frequency, frequency_levels, [1.] + parameters_common[6:9])
        par[6] -= np.interp(frequency, frequency_levels, [0] + parameters_common[14:17])
        return par


def get_parameters_for_frequency(f):
    if f == "common":
        return parameters_common
    else:
        par = np.zeros(19)
        par[6:9] = 1
        par[0:6] = parameters[f][0:6]
        par[9:14] = parameters[f][6:11]
        par[17:19] = parameters[f][11:13]
        return par

# parameters_common = [2.20517242e+00, 3.72924098e-08, 3.77027776e-02, 7.16206717e-04,
#                     1.34590082e+01, 2.12620654e+02, 8.65661146e-01, 5.17534627e-01,
#                     3.43653164e-01,   8.67280709e-02,  3.86860955e+01, -8.28949637e+04,
#                      -3.60586924e+00, -5.46979329e+01, -1.14237778e-01, -6.45139549e-01,
#                      7.93158176e-01,  5.49430755e-03, -7.83909469e-05,  1.17499633e-03]

#  9.54421382e-02  2.69354907e+01 -8.01765492e+03 -2.92460195e+00
#  -3.37589449e+00  1.80937949e-02 -1.43360431e-01  1.28693865e-01
#   4.88868664e-03 -6.12340751e-05  5.87790419e-03


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
