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

x = np.linspace(0, 1, 1000)
z_vec = np.linspace(0, 1, 10)

for z in z_vec:
    r = np.sqrt(x**2 + z**2)
    s = z**3/r**5
    print(s)
    plt.plot(x, s)

plt.show()