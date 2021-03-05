import os
import pickle
import subprocess

import numpy as np

from common import abq_viewer


def get_data_from_path(path_points, odb_filename, variable, component=None, step_name=None, frame_number=None,
                       output_position='ELEMENT_NODAL'):
    work_directory = os.path.dirname(os.path.abspath(odb_filename))
    parameter_pickle_name = work_directory + '/parameter_pickle.pkl'
    path_points_filename = work_directory + '/path_points.npy'
    data_filename = work_directory + '/path_data.npy'
    parameter_dict = {'odb_filename': odb_filename,
                      'variable': variable,
                      'path_points_filename': path_points_filename,
                      'data_filename': data_filename}
    if component is not None:
        parameter_dict['component'] = component
    if step_name is not None:
        parameter_dict['step_name'] = step_name
    if frame_number is not None:
        parameter_dict['frame_number'] = frame_number
    parameter_dict['output_position'] = output_position

    with open(parameter_pickle_name, 'wb') as pickle_file:
        pickle.dump(parameter_dict, pickle_file, protocol=2)
    if not isinstance(path_points, np.ndarray):
        path_points = np.array(path_points)
    np.save(path_points_filename, path_points)
    os.chdir('abaqus_functions')
    print(abq_viewer + ' noGUI=write_data_along_path.py -- ' + parameter_pickle_name)
    job = subprocess.Popen(abq_viewer + ' noGUI=write_data_along_path.py -- ' + parameter_pickle_name,
                           shell=True)
    dfgfdgdfgdg
    job.wait()
    os.chdir('..')
    data = np.unique(np.load(data_filename), axis=0)
    os.remove(parameter_pickle_name)
    os.remove(path_points_filename)
    os.remove(data_filename)
    _, idx = np.unique(data[:, 0], return_index=True)
    return data[idx, 1]


def main():
    odb_file_name = os.path.expanduser('~/railway_ballast/python/embankment_model/results_2.odb')
    path_points = np.zeros((100, 3))
    path_points[:, 0] = np.linspace(0.1, 0.2, 100)
    path_points[:, 1] += 8.
    path_points[:, 2] += 1.

    print(get_data_from_path(path_points, odb_file_name, 'UP', 'UP2', output_position='NODAL'))


if __name__ == '__main__':
    main()
