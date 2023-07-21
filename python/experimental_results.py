import glob
import os
import pathlib
import re

from numbers import Number

import numpy as np


class Experiment:
    def __init__(self, p, q, f, filename_axial, filename_volumetric=None, num_points=100):
        self.p = p
        self.q = q
        self.f = f

        # Reading axial data
        data = np.genfromtxt(filename_axial, delimiter=',')
        idx = np.argsort(data[:, 0])
        data = data[idx, :]
        axial_cycles = np.log(data[:, 0])
        axial_strain = data[:, 1]

        if filename_volumetric is None:
            filename_volumetric = 'volumetric_strain_' + filename_axial.lstrip('axial_strain_')
        data = np.genfromtxt(filename_volumetric, delimiter=',')
        idx = np.argsort(data[:, 0])
        data = data[idx, :]
        volumetric_cycles = np.log(data[:, 0])
        volumetric_strain = data[:, 1]
        self.cycles = np.linspace(np.max([axial_cycles[0], volumetric_cycles[0]]),
                                  np.min([axial_cycles[-1], volumetric_cycles[-1]]), num_points)

        self.axial_strain = 0*self.cycles + 1 + axial_strain[0]
        self.volumetric_strain = 0*self.cycles + 1. + volumetric_strain[0]
        self.axial_strain = np.interp(self.cycles, axial_cycles, axial_strain)
        self.volumetric_strain = np.interp(self.cycles, volumetric_cycles, volumetric_strain)
        self.cycles = np.exp(self.cycles)

    def deviatoric_axial_strain(self):
        edev = self.axial_strain - self.volumetric_strain/3
        return edev


class ExperimentalResults:
    def __init__(self):
        self.data = []

    def read(self, directory, regex_str=r'axial_strain_p=(\d+)kPa_q=(\d+)kPa_f=(\d+)Hz.dat'):
        curr_dir = os.getcwd()
        os.chdir(directory)
        filenames = glob.glob('*.dat')
        for filename in filenames:
            match = re.search(regex_str, filename)
            if match:
                self.data.append(Experiment(p=float(match.group(1)), q=float(match.group(2)), f=float(match.group(3)),
                                            filename_axial=match.group(0)))
        os.chdir(curr_dir)

    def get_data(self, p=None, q=None, f=None):
        if p is None and q is None and f is None:
            return self.data
        sub_set = []
        if isinstance(p, Number):
            p = [p]
        if isinstance(q, Number):
            q = [q]
        if isinstance(f, Number):
            f = [f]
        for dataset in self.data:
            if p is None or dataset.p in p:
                if q is None or dataset.q in q:
                    if f is None or dataset.f in f:
                        sub_set.append(dataset)
        return sub_set


sun_et_al_16 = ExperimentalResults()
sun_et_al_16.read(os.path.expanduser('~/railway_ballast/experimental_data/sun_et_al_16/'))

fouled_folder = pathlib.Path('~/railway_ballast/experimental_data/tennakoon_indraratna_14/').expanduser()
fouled = {
    0: Experiment(10, 230, 20, fouled_folder / 'axial_0_vci.csv', fouled_folder / 'vol_0_vci.csv'),
    10: Experiment(10, 230, 20, fouled_folder / 'axial_10_vci.csv', fouled_folder / 'vol_10_vci.csv'),
    25: Experiment(10, 230, 20, fouled_folder / 'axial_25_vci.csv', fouled_folder / 'vol_25_vci.csv'),
    50: Experiment(10, 230, 20, fouled_folder / 'axial_50_vci.csv', fouled_folder / 'vol_50_vci.csv'),
    80: Experiment(10, 230, 20, fouled_folder / 'axial_80_vci.csv', fouled_folder / 'vol_80_vci.csv')
}

def main():
    from collections import namedtuple
    import matplotlib.pyplot as plt
    import matplotlib.style
    matplotlib.style.use('classic')
    plt.rc('text', usetex=True)
    plt.rc('font', serif='Computer Modern Roman')
    plt.rcParams.update({'font.size': 20})
    plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
    plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                      'monospace': ['Computer Modern Typewriter']})

    symbols = {5: 's', 10: 'o', 20: '^', 30: 'v', 40: 'D', 50: '<', 60: '>'}

    Figure = namedtuple('Figure', ['number', 'p', 'q'])
    figures = [Figure(number=0, p=[10], q=[230]),
               Figure(number=1, p=[30], q=[230, 276]),
               Figure(number=2, p=[60], q=[230]),
               Figure(number=3, p=[60], q=[370., 460.])]
    for figure in figures:

        experimental_data = sun_et_al_16.get_data(p=figure.p, q=figure.q)
        for data_set in experimental_data:
            plt.figure(2*figure.number)
            plt.semilogx(data_set.cycles, data_set.axial_strain, '-' + symbols[data_set.f], lw=2)
            plt.ylim(0, 0.35)
            plt.figure(2*figure.number+1)
            plt.semilogx(data_set.cycles, data_set.volumetric_strain, '-' + symbols[data_set.f], lw=2)
            plt.ylim(-0.06, 0.24)
    plt.show()


if __name__ == '__main__':
    main()
