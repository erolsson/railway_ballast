import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style

from finite_element_model.simulations import simulation1
from get_data_from_path import get_data_from_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})


def main():
    results_odb_filename = os.path.expanduser('//results_2.odb')
    ballast_start_height = 0
    for layer in simulation1.layers:
        if layer.name.lower().startswith('ballast'):
            break
        ballast_start_height += layer.height
    total_height = sum([layer.height for layer in simulation1.layers])
    path_points = np.zeros((1000, 3))
    y = np.linspace(total_height, ballast_start_height, 1000)
    path_points[:, 1] = y
    path_points[:, 0] = simulation1.track_gauge/2
    path_points[:, 2] += 1e-6
    for n in simulation1.cycles:
        step_name = 'cycles_' + str(int(n))
        up = get_data_from_path(path_points, results_odb_filename, 'UP', 'UP2', output_position='NODAL',
                                step_name=step_name)
        plt.plot(y[0] - y, -up*1000, lw=2, label='$N=' + str(int(n)) + '$')
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Settlement [mm]', fontsize=24)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('settlement.png')
    plt.show()


if __name__ == '__main__':
    main()
