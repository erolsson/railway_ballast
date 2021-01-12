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

parameters[5] = [1.18705105e+00,  3.75814744e-04, 5.43522967e+00,  4.54045581e-02,
                 2.94230799e+01,  8.08711944e+00]

parameters[10] = [6.84064006e-01,  2.87438142e-04, 5.43522967e+00,  4.54045581e-02,
                  2.49613427e+01,  8.07865840e+00]

parameters[20] = [9.39141439e-01,  1.90573525e-04, 4.91382134e+00,  6.19133427e-02,
                  1.09092425e+01,  8.01956620e+00]

parameters[40] = [1.06807362e+00, 2.03734495e-04, 5.43068531e+00, 4.08554859e-02,
                  6.96762283e+00, 8.93961830e+00]

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
