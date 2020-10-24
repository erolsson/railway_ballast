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

par = {5: np.array([4.29466804e+00, 7.95338549e-13, -2.49079302e+01, -2.63383075e-03,
                    1.40385907e+01, 2.67798569e+02, 9.71329457e-01, 1.41773676e+00,
                    2.07605043e+00]),

       10: np.array([3.25230313e+00, 1.74540229e-13, -7.67149741e+01, -1.13863239e-02,
                     1.35060633e+01, 1.73030452e+03, 1.42144129e-01, 2.16274731e-01,
                     3.02304057e-01]),

       20: np.array([1.92510062e+00, 1.26427910e-11, 6.67370258e+03, -1.44282399e-01,
                     2.75744815e+00, 3.32909517e+04, 1.60779590e-02, 2.49632584e-02,
                     3.47000750e-02]),

       40: np.array([1.92349523e+00, 2.97337921e-11, - 4.73646242e+04, - 1.03805831e+00,
                     1.34217785e+00, 3.22079739e+04, 2.42366633e-02, 3.22609229e-02,
                     4.55613347e-02])}

pressures = np.array([10, 30, 60])
ep = np.linspace(0, 0.3, 1000)
for f, parameters in par.items():
    plt.figure(0)
    p_par = parameters[6:9]
    plt.plot(pressures, p_par/p_par[0])

    plt.figure(1)
    plt.plot(ep, parameters[5]*(1-np.exp(-parameters[4]*ep)))
    plt.plot()

plt.figure(0)
plt.plot(pressures, 0.58*np.sqrt(1+0.2*pressures+0*pressures**2), '--')
plt.show()
