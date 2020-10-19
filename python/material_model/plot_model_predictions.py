import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from experimental_results.experimental_results import sun_et_al_16
from functional_shapes import permanent_strain

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


# f = 5 Hz

par = np.array([3.18452468e+00,  5.79153855e-10,  8.13999845e-01, -1.57228447e-04,
                1.68891350e+01, -2.00022511e+01,  1.05777046e+02])

# par = np.array([3.13428661e+00,  3.62814163e-10,  4.52625876e-01,  4.26854206e-10,
#                1.75132361e+01, -2.18733820e+01,  1.30768783e+02])

# par = np.array([3.34868394e+00,  4.65718114e-10,  1.00535319e+00,  1.07864924e+00,
#                1.38303647e-10,  1.62933421e+01, -1.72924629e+01,  7.61020870e+01])

# f = 10 Hz

par = np.array([1.81886025e+00,  8.09590763e-09,  3.02248576e-02,  6.16678475e-04,
                1.94384809e+01, -1.72760368e+00,  5.06221378e+00])
# par = np.array([3.65260909e+00,  2.12789850e-11,  1.96867169e-01,  8.06196338e-10,
#                1.34770968e+01, -2.51841032e+01,  1.89992871e+02])

# par = np.array([9.57058472e-01,  7.49776116e-07,  1.27648470e-01,  4.74085636e-05,
#                2.85328510e+01, -1.01461558e+02,  2.55740255e+02])

# par = np.array([3.65260891e+00,  6.78995846e-11,  1.88764564e+00,  3.71616109e-01,
#                1.77354432e-10,  1.34771050e+01, -1.83301606e+01,  1.19955349e+02])


# f = 20 Hz
par = np.array([ 1.35966478e+00,  2.68774845e-05,  2.25879758e+03, -3.81948760e-03,
                 1.08065068e+00, -9.13254681e-02,  1.42448926e+01,  1.26502963e+00,
                 1.80823494e+00,  2.45840273e+00])

par = np.array([ 1.48793286e+00,  2.23887423e-08,  3.47498051e-02, -1.31835549e-03,
                1.05564935e+00, -1.30507801e-01,  1.52408961e+01,   1.24832180e+00,  1.76638543e+00,  2.29595102e+00])

# f = 40 Hz
# par = np.array([1.90166096e+00,  1.39154106e-09,  1.86960565e-01,  2.69030714e-02,
#                -1.48790476e-10,  3.42180958e+00, -1.31346021e+02,  6.47658273e+02])

par = np.array(par)
for i, freq in enumerate([5., 10., 20., 40]):
    experimental_data = sun_et_al_16.get_data(f=freq)
    for experiment in experimental_data:
        plt.figure(i)
        plt.semilogx(experiment.cycles, experiment.strain)
        e0 = experiment.strain[0]
        cycles = np.exp(np.linspace(np.log(1), np.log(5e5), 100))
        model_strain = permanent_strain(cycles, experiment.p, experiment.q, freq, par, e0=0)
        plt.semilogx(cycles, model_strain + e0, '--x')
    plt.figure(i)
    plt.xlim(1, 5e5)
    plt.ylim(0., 0.3)

plt.show()
