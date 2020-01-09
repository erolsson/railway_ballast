import imp
import os
import pickle
import subprocess


with open('python_modules.pkl') as module_pickle:
    modules = pickle.load(module_pickle)

numpy = imp.load_module('numpy', modules['numpy'][0], modules['numpy'][1], modules['numpy'][2])
scipy = imp.load_module('scipy', modules['scipy'][0], modules['scipy'][1], modules['scipy'][2])


