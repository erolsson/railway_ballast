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
                 1.25454110e+01, 1.94706613e+02, 1.46936052e-01,  4.51955283e+01,
                 -5.12030540e+03, -7.15689757e+00, -1.27848215e+00,  8.65098143e-03,
                 -1.06993520e-04]

parameters[10] = [1.81696447e+00, 3.67922225e-07, 5.25068390e-02, 7.49368438e-04,
                  1.39169448e+01, 1.86223864e+02,  7.73461328e-02,  4.71392820e+01,
                  -9.82332826e+03, -6.38557798e+00, -5.04045480e+00,  6.22446188e-03,
                  -8.26719358e-05]

parameters[20] = [1.41053055e+00, 3.32379026e-06, 8.67794008e-01, 4.52482073e-04,
                  2.00671901e+00, 1.65727190e+02, 6.80413393e-02, -1.14943221e+01,
                  5.75442841e+03, -4.14366403e-02, -1.49249306e+00, 6.74729195e-03,
                  -7.82340337e-05]

parameters[40] = [1.40956811e+00, 1.10662626e-06, 7.38373606e-02, 4.52749231e-04,
                  1.25748667e+00, 4.76597957e+02, 3.33608109e-03, -1.20218444e+01,
                  2.98562054e+03, -1.04935998e-01, -1.80513397e+00,  2.62758004e-03,
                  -3.03834857e-05]

parameters_common = [2.20517242e+00, 3.72924098e-08, 3.77027776e-02, 7.16206717e-04,
                    1.34590082e+01, 2.12620654e+02, 8.65661146e-01, 5.17534627e-01,
                    3.43653164e-01,   1.35237246e-01,  2.11811921e+01, -7.13402604e+03,
                     -2.00586879e+00, -2.82053712e+00,  3.81672072e-02, -1.03244401e-01,
                     1.32865400e-01,  4.18561848e-03, -4.99218275e-05]
# 3.60519179e-06  2.42571225e-06  4.69929518e+03  4.02427403e-01
#   2.16544726e-02  1.00335940e+00  2.99206928e+00  6.26584857e+00
#   1.72037575e-03 -2.50345880e-03]


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
