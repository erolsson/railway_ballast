import numpy as np

from material_model import MaterialModel
from model_parameters import get_parameters

static_stress = np.array([2, -32.5, -12, 0, 0, 0])
cyclic_stress = np.array([6, -5, 5, 0, 0, 0])
parameters = get_parameters(40, common=False)
model = MaterialModel(parameters)
model.update(np.array([0, 1e6]), cyclic_stress, static_stress)
print(model.strain()[-1, :])