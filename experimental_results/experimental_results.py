import glob
import os
import re

from numbers import Number

import numpy as np


class Experiment:
    def __init__(self, p, q, f, filename):
        self.p = p
        self.q = q
        self.f = f

        data = np.genfromtxt(filename, delimiter=',')
        if data[-1, 0] < data[0, 0]:
            data = np.flipud(data)
        self.cycles = data[:, 0]
        self.strain = data[:, 1]


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
sun_et_al_16.read(os.path.expanduser('~/work/railway_ballast/experimental_results/sun_et_al_16'))
