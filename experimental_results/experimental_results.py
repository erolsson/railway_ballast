import glob
import os
import re

from numbers import Number

import numpy as np


class Experiment:
    def __init__(self, p, q, f, filename, num_points=100):
        self.p = p
        self.q = q
        self.f = f

        # Reading axial data
        data = np.genfromtxt(filename, delimiter=',')
        idx = np.argsort(data[:, 0])
        data = data[idx, :]
        axial_cycles = np.log(data[:, 0])
        axial_strain = data[:, 1]

        volumetric_filename = 'volumetric_strain_' + filename.lstrip('axial_strain_')
        data = np.genfromtxt(volumetric_filename, delimiter=',')
        idx = np.argsort(data[:, 0])
        data = data[idx, :]
        volumetric_cycles = np.log(data[:, 0])
        volumetric_strain = data[:, 1]

        self.cycles = np.linspace(np.max([axial_cycles[0], volumetric_cycles[0]]),
                                  np.min([axial_cycles[-1], volumetric_cycles[-1]]),
                                  num_points)

        self.axial_strain = np.interp(self.cycles, axial_cycles, axial_strain)
        self.volumetric_strain = np.interp(self.cycles, volumetric_cycles, volumetric_strain)
        self.cycles = np.exp(self.cycles)

    def deviatoric_axial_strain(self):
        return self.axial_strain - self.volumetric_strain/3


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
                                            filename=match.group(0)))
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
sun_et_al_16.read(os.path.expanduser('~/railway_ballast/experimental_results/sun_et_al_16'))
