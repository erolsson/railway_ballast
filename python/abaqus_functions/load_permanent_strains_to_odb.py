import sys
import pickle

import numpy as np

from odb_io_functions import write_field_to_odb


def write_permanent_strains_to_odb(odb_file, array_file_name, pickle_file_name):
    ep = np.load(array_file_name)
    with open(pickle_file_name, 'rb') as pickle_file:
        data = pickle.load(pickle_file)
    instance_name = data['instance']
    print(instance_name)
    element_set_name = data['element_set']
    cycles = data['cycles']
    for i, n in enumerate(cycles):
        step_name = 'cycles=' + str(int(n))
        print(ep[i, :, :].shape)
        write_field_to_odb(ep[i, :, :], 'EP', odb_file, step_name=step_name, instance_name=instance_name,
                           set_name=element_set_name, step_description='Results after ' + str(n) + ' cycles',
                           field_description='Permanent strain')


if __name__ == '__main__':
    write_permanent_strains_to_odb(sys.argv[-3], sys.argv[-2], sys.argv[-1])
