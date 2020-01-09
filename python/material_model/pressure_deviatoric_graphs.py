import matplotlib.pyplot as plt
import matplotlib.style

import numpy as np

from functional_shapes import permanent_strain

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}", r"\usepackage{bm}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

p = 30
q_vals = np.arange(100, 600, 100)
n = np.exp(np.linspace(0, np.log(1e6)))

par = np.array([5.19502802e+00, 7.64481869e-07, 4.70091315e+01, 1.34948445e-13,
                1.14929617e+01, 1.12439205e+01])

plt.figure(0)
for q in q_vals:
    e = permanent_strain(n, p, q, par)
    plt.semilogx(n, e, '--', lw=2)

plt.text(1000, 0.23, r'{$\bm p=$ \textbf{30 kPa}}')

plt.figure(1)
q = 300
p_vals = np.arange(10, 60, 10)
for p in p_vals:
    e = permanent_strain(n, p, q, par)
    plt.semilogx(n, e, '--', lw=2)

plt.text(1000, 0.23, r'{$\bm q=$ \textbf{300 kPa}}')

for fig in [0, 1]:
    plt.figure(fig)
    plt.ylim(0, 0.25)
    plt.xlabel('Cycles $N$')
    plt.ylabel(r'Permanent strain $\varepsilon_p$')
    plt.tight_layout()
    
plt.figure(0)
plt.savefig('../Figures/pressure.png')

plt.figure(1)
plt.savefig('../Figures/deviatoric.png')
plt.show()
