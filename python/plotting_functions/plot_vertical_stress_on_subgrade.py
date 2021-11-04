import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from plotting_functions.plot_stresses_along_path import get_data_from_path
from comparison_of_models import get_path_points_for_fem_simulation


matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
figure_directory = os.path.expanduser('~/railway_ballast/Figures/')
load_levels = np.array([17.5, 22.5, 30.0])

for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
    for geometry, color in zip(['low', 'high'], ['r', 'b']):
        vertical_stress = 0*load_levels
        for i, load in enumerate(load_levels):
            path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
            odb_filename = (odb_directory + '/stresses_' + rail_fixture + '_' + geometry + '_'
                            + str(load).replace('.', '_') + 't.odb')
            s = (get_data_from_path(path_points, odb_filename, 'S', step_name='gravity', component='S22',
                                    output_position='INTEGRATION_POINT') +
                 get_data_from_path(path_points, odb_filename, 'S', step_name='cyclic_stresses', component='S22',
                                    output_position='INTEGRATION_POINT'))
            vertical_stress[i] = s[-1]
        plt.plot(load_levels, -vertical_stress/1000, line + color, lw=2,
                 label=geometry.capitalize() + ' embankment, ' + rail_fixture)

plt.xlabel('Axle load [tonnes]', fontsize=24)
plt.ylabel('Compressive stress on subgrade [kPa]', fontsize=24)
plt.legend(loc='best')
plt.xlim(17.5, 30.)
plt.tight_layout()
plt.savefig(figure_directory + '/subgrade_stresses.png')
plt.show()
