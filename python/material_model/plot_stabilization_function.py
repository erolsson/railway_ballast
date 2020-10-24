import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import fmin

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def hf(ep):
    # nf = parameters[4]
    # H0 = parameters[5]
    # H1 = parameters[6]
    # return H0 + H1*(1-np.exp(-nf*(ep - e0)))
    f_vals = [8.13220695e-02,
              2.68868736e+01, 4.91914262e+01, -2.45166693e+01, -2.33720631e+01,
              -4.08586276e+01, -1.00691061e-02, -4.69099710e-01, 6.33116729e+01,
              -4.92231884e-01, -5.99374859e+01, 4.48312718e+00, -4.41065771e-02,
              -5.29131981e+00, -1.17381751e+02, 4.53602062e+00, -9.32487550e+00,
              9.18212675e-02, 8.05041126e+00, 1.19068324e+00]

    strain_values = np.linspace(0, 0.3, 20)
    f_vals = np.abs(f_vals)
    return np.interp(ep, strain_values, np.cumsum(f_vals))


def model(par, ep):
    return par[0]*(1-np.exp(-par[1]*ep))


def residual(par, ep, y):
    return np.sum((model(par, ep) - y)**2)


def main():
    ep1 = np.linspace(0, 0.3, 1000)
    plt.plot(ep1, 4.47862581e+04*(1-np.exp(-2.91700112e-02*ep1)))
    plt.plot(ep1, 350*(1 - np.exp(-10*ep1)))
    plt.show()


if __name__ == '__main__':
    main()
