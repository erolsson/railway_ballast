import glob
import re

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from scipy.optimize import curve_fit

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

groups = [{'p': 10, 'q': [230]},
          {'p': 30, 'q': [230, 276]},
          {'p': 60, 'q': [230]},
          {'p': 60, 'q': [370, 460]}]

frequencies_symbols = {5: 'bs', 10: 'go', 20: 'r^', 30: 'kv', 40: 'md', 50: 'y<', 60: 'b>'}

frequency_data = []
for i, test_series in enumerate(groups):
    for q in test_series['q']:
        p = int(test_series['p'])
        filename_template = 'sun_et_al_16/axial_strain_p=' + str(test_series['p']) + 'kPa_q=' + str(q) + \
                            'kPa_f=*Hz.dat'
        filenames = glob.glob(filename_template)

        for filename in filenames:
            f = re.match(filename_template.replace('*', r'(\d+)'), filename).groups()[0]
            if p == 60:
                filename = filename.replace('*', f)
                data = np.genfromtxt(filename, delimiter=',')
                if data[0, 0] > data[-1, 0]:
                    data = np.flipud(data)
                plt.figure(1)

                if data[-1, 0] > 2e5 and np.diff(data[:, 1])[-1] < 1e-3:
                    plt.plot(float(f), data[-1, 1], 'x', mew=2)
                    plt.figure(2)
                    if q == 230.:
                        plt.plot(float(f), data[-1, 1]/0.04, 'x', mew=2)
                        frequency_data.append([float(f), data[-1, 1]/0.04])
                    if q == 370.:
                        plt.plot(float(f), data[-1, 1]/0.077, 'x', mew=2)
                        frequency_data.append([float(f), data[-1, 1]/0.077])
                e1000 = np.interp(1000, data[:, 0], data[:, 1])

frequency_data = np.array(frequency_data)


def power(x, a, b, c):
    return 1 + a*(1 - np.exp(-b*x**c))


par = curve_fit(power, frequency_data[:, 0], frequency_data[:, 1], maxfev=int(1e8))
f = np.linspace(0., 100., 1000)
print par
plt.figure(2)
plt.plot(f, power(f, par[0][0], par[0][1], par[0][2]))
plt.show()
