import os
import subprocess

from common import abq
from read_data_from_odb import read_data_from_odb
from write_data_to_odb import write_data_to_odb


odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
ballast_element_set = 'EMBANKMENT_INSTANCE_BALLAST_ELEMENTS'


def main():
    loads = [17.5]
    for geometry in ['low', 'high']:
        for rail_fixture in ['sleepers', 'slab']:
            for load in loads:
                stress_state_odb_filename = (odb_directory + '/stresses_' + rail_fixture + '_' + geometry + '_'
                                             + str(load).replace('.', '_') + 't.odb')
                sim_odb_filename = (odb_directory + '/embankment_' + rail_fixture + '_' + geometry + '_'
                                    + str(load).replace('.', '_') + 't.odb')
                static_stresses = read_data_from_odb('S', sim_odb_filename, step_name='gravity',
                                                     set_name=ballast_element_set)
                if os.path.isfile(stress_state_odb_filename):
                    os.remove(stress_state_odb_filename)
                os.chdir('abaqus_functions')
                job = subprocess.Popen(abq + ' python create_empty_odb.py ' + stress_state_odb_filename + ' ' +
                                       sim_odb_filename,  shell=True)
                job.wait()
                os.chdir('..')

                write_data_to_odb(static_stresses, 'S', stress_state_odb_filename, set_name=ballast_element_set,
                                  step_name='gravity')
                max_stresses = read_data_from_odb('S', sim_odb_filename, step_name='loading',
                                                  set_name=ballast_element_set)
                cyclic_stresses = max_stresses - static_stresses
                write_data_to_odb(cyclic_stresses, 'S', stress_state_odb_filename, set_name=ballast_element_set,
                                  step_name='cyclic_stresses')


if __name__ == '__main__':
    main()
