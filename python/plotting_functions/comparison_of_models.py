import numpy as np
from finite_element_model.simulations import simulations
def get_path_points_for_fem_simulation(sim_name):
    ballast_start_height = 0
    fem_simulation = simulations[sim_name]
    for layer in fem_simulation.layers:
        if layer.name.lower().startswith('ballast'):
            break
        ballast_start_height += layer.height
    total_height = sum([layer.height for layer in fem_simulation.layers])

    path_points = np.zeros((1000, 3))
    y = np.linspace(total_height-1e-3, ballast_start_height+1e-3, 1000)
    path_points[:, 1] = y
    path_points[:, 0] = fem_simulation.track_gauge/2
    path_points[:, 2] += 1e-3

    return path_points