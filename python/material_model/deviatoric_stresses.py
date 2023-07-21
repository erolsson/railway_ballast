import numpy as np

import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.gridspec as gridspec

from multiprocesser.multiprocesser import multi_processer

from experimental_results import fouled
from material_model import MaterialModel
from model_parameters import get_parameters


def main():
    frequencies = [5, 10, 20, 40]
    pressures = [30, 60]
    colors = ['r', 'b', 'g', 'k']
    deviatoric_stresses = np.linspace(1, 500, 100)
    n = np.exp(np.linspace(0, np.log(5e5)))

    for p in pressures:
        for c, f in zip(colors, frequencies):
            static_stress = -p*np.array([1, 1, 1, 0, 0, 0])
            strain = 0*deviatoric_stresses
            for i, q in enumerate(deviatoric_stresses):
                cyclic_stress = -q*np.array([1, 0, 0, 0, 0, 0])
                parameters = get_parameters(f, common=False)
                model = MaterialModel(parameters)
                model.update(n, cyclic_stress, static_stress)
                strain[i] = -model.deviatoric_strain()[-1, 0]
            plt.plot(deviatoric_stresses, strain, c, lw=2)

    plt.ylim(0, 0.3)
    plt.show()


if __name__ == '__main__':
    main()
