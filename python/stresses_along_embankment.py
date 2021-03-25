import os

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.style

from get_data_from_path import get_data_from_path
from comparison_of_models import get_path_points_for_fem_simulation
from plot_stresses_along_path import mises, get_tensor_from_path

matplotlib.style.use('classic')
plt.rc('text', usetex=True)
plt.rc('font', serif='Computer Modern Roman')
plt.rcParams.update({'font.size': 20})
plt.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'],
                  'monospace': ['Computer Modern Typewriter']})

odb_directory = os.path.expanduser('~/railway_ballast/odbs/')


def main():
    plt.figure(0, figsize=(12, 9))
    gs = gridspec.GridSpec(3, 2)
    ax1 = plt.subplot(gs[0:2, 0:1])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Cyclic von Mises stress, $q$ [kPa]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(0, 30)
    ax1.yaxis.set_label_coords(-0.12, 0.5)
    plt.tight_layout()

    ax2 = plt.subplot(gs[0:2, 1:2])
    plt.xlabel('Distance from ballast surface [m]', fontsize=24)
    plt.ylabel('Cyclic von Mises stress, $q$ [kPa]', fontsize=24)
    plt.xlim(0, 4.3)
    plt.ylim(0, 150)
    ax2.yaxis.set_label_coords(-0.12, 0.5)
    plt.tight_layout()

    for i, rail_fixture, line in enumerate(['slab', 'sleepers']):
        path_points = get_path_points_for_fem_simulation(rail_fixture + '_high')
        for load, c in zip([17.5, 22.5, 30.], ['g', 'r', 'b']):
            pass

