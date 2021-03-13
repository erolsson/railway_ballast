import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from finite_element_model.simulations import simulations
from get_data_from_path import get_data_from_path
from comparison_of_models import get_path_points_for_fem_simulation

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/python/embankment_model')

for geometry in ['low', 'high']:
    for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
        for load, c in zip([22.5, 30.], ['r', 'b']):
            odb_file_name = (odb_directory + '/embankment_' + rail_fixture + '_' + geometry + '_'
                             + str(load).replace('.', '_') + 't.odb')
            print(odb_file_name)
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            stress = get_data_from_path(path_points, odb_file_name, 'S', component='S11')
            print(stress)