from __future__ import print_function

import imp
import pickle
import sys

try:
    with open('python_modules.pkl') as module_pickle:
        modules = pickle.load(module_pickle)
except IOError:
    print("python_modules not set!")
    print("Run \"python find_modules.py\" from the terminal")
    sys.exit()

numpy = imp.load_module('numpy', modules['numpy'][0], modules['numpy'][1], modules['numpy'][2])
scipy = imp.load_module('scipy', modules['scipy'][0], modules['scipy'][1], modules['scipy'][2])


