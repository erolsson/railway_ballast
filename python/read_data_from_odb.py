import os
import pickle
import subprocess

from common import abq, create_temp_dir_name


def read_data_from_odb(field_id, odb_file_name, step_name=None, frame_number=-1, set_name='', instance_name='',
                       get_position_numbers=False, get_frame_value=False, position='INTEGRATION_POINT'):

    work_directory = create_temp_dir_name(odb_file_name)
    os.makedirs(work_directory)
    parameter_pickle_name = work_directory + '/parameter_pickle.pkl'
    results_pickle_name = work_directory + '/results.pkl'
    with open(parameter_pickle_name, 'wb') as pickle_file:
        pickle.dump({'field_id': field_id, 'odb_file_name': odb_file_name, 'step_name': step_name,
                     'frame_number': frame_number, 'set_name': set_name, 'instance_name': instance_name,
                     'get_position_numbers': get_position_numbers, 'get_frame_value': get_frame_value,
                     'position': position},
                    pickle_file, protocol=2)
    os.chdir('abaqus_functions')
    job = subprocess.Popen(abq + ' python read_data_from_odb.py ' + parameter_pickle_name + ' ' + results_pickle_name,
                           shell=True)
    job.wait()
    os.chdir('..')
    with open(results_pickle_name, 'rb') as results_pickle:
        data = pickle.load(results_pickle, encoding='latin1')
    os.remove(parameter_pickle_name)
    os.remove(results_pickle_name)
    os.removedirs(work_directory)
    if not get_position_numbers and not get_frame_value:
        return data['data']
    elif not get_position_numbers:
        return data['data'], data['frame_value']
    elif not get_frame_value:
        return data['data'], data['node_labels'], data['element_labels']
    else:
        return data['data'], data['frame_value'], data['node_labels'], data['element_labels']
