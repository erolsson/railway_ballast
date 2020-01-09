import imp
import os
import pickle
import subprocess


work_dir = os.getcwd()
os.chdir(os.path.abspath(os.path.dirname(__file__)))
subprocess.Popen('python find_numpy.py', shell=True)
with open('python_modules.pkl') as module_pickle:
    modules = pickle.load(module_pickle)
print modules
numpy = imp.load_module('numpy', modules['numpy'][0], modules['numpy'][1], modules['numpy'][2])
scipy = imp.load_module('scipy', modules['scipy'][0], modules['scipy'][1], modules['scipy'][2])
print numpy.__version__, scipy.__version__

os.chdir(work_dir)
