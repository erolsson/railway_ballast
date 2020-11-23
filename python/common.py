from __future__ import print_function
import os
import imp
import pickle
import sys

package_path = os.path.dirname(__file__)

try:
    with open(package_path + '/python_modules.pkl') as module_pickle:
        modules = pickle.load(module_pickle)
except IOError:
    print("python_modules not set!")
    print("Run \"python find_modules.py\" from the terminal")
    sys.exit()

sys.path.append()
six = imp.load_module('six', open(modules['six'][0]), modules['six'][1], modules['six'][2])
numpy = imp.load_module('numpy', modules['numpy'][0], modules['numpy'][1], modules['numpy'][2])
scipy = imp.load_module('scipy', modules['scipy'][0], modules['scipy'][1], modules['scipy'][2])
integrate = imp.load_module('scipy.integrate', modules['scipy'][0], modules['scipy'][1], modules['scipy'][2])