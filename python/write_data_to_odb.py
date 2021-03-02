import os
import pickle
import subprocess

import numpy as np

from common import abq


def write_data_to_odb(field_data, field_id, odb_file_name, step_name, instance_name='', set_name='',
                      step_description='', frame_number=None, frame_value=None, field_description='',
                      position='INTEGRATION_POINT'):
    work_directory = os.path.splitext(odb_file_name)[0] + '_tempdir'
    os.makedirs(work_directory)
    pickle_filename = work_directory + r'/load_field_to_odb_pickle.pkl'
    data_filename = work_directory + r'/field_data.npy'
    np.save(data_filename, field_data)

    with open(pickle_filename, 'wb') as pickle_file:
        pickle.dump({'field_id': field_id, 'odb_file': odb_file_name, 'step_name': step_name,
                     'instance_name': instance_name, 'set_name': set_name, 'step_description': step_description,
                     'frame_number': frame_number, 'frame_value': frame_value, 'field_description': field_description,
                     'position': position},
                    pickle_file, protocol=2)

    os.chdir('abaqus_functions')
    job = subprocess.Popen(abq + ' python write_data_to_odb.py ' + data_filename + ' ' + pickle_filename, shell=True)
    job.wait()
    os.chdir('..')

    os.remove(pickle_filename)
    os.remove(data_filename)
    os.removedirs(work_directory)
