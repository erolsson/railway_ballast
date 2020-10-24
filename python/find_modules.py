import imp
import pickle

python_modules = {name: imp.find_module(name) for name in ['numpy', 'scipy']}
python_modules['six'] = (imp.find_module('six')[-2], imp.find_module('six')[-2], imp.find_module('six')[-1])
with open('python_modules.pkl', 'w') as pickle_handle:
    pickle.dump(python_modules, pickle_handle)
