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


def get_stress_tensor_from_path(odb_file_name, path_points, step_name=None, frame_number=None):
    stress_components = ['S11', 'S22', 'S33', 'S12', 'S13', 'S23']
    data = np.zeros((path_points.shape[0], 6))
    for i, component in enumerate(stress_components):
        stress = get_data_from_path(path_points, odb_file_name, 'S', step_name=step_name,
                                    frame_number=frame_number, output_position='INTEGRATION_POINT',
                                    component=component)
        data[:, i] = stress
    return data


def main():
    for geometry in ['low', 'high']:
        for rail_fixture, line in zip(['slab', 'sleepers'], ['--', '-']):
            min_stresses = None
            for load, c in zip([22.5, 30.], ['r', 'b']):
                odb_file_name = (odb_directory + '/embankment_' + rail_fixture + '_' + geometry + '_'
                                 + str(load).replace('.', '_') + 't.odb')

                path_points = get_path_points_for_fem_simulation(rail_fixture + '_' + geometry)
                if load == 22.5:
                    min_stresses = get_stress_tensor_from_path(odb_file_name, path_points, step_name='gravity')
                    static_pressure = -np.sum(min_stresses[:, :3], axis=1)/3
                    plt.plot(path_points[0, 1] - path_points[:, 1], static_pressure, line + c, lw=2)
                # max_stresses = get_stress_tensor_from_path(odb_file_name, path_points)
    plt.show()


if __name__ == '__main__':
    main()