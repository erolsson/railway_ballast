import odbAccess

import pickle
import sys

from odb_io_functions import read_field_from_odb


def write_stress_pickles(stress_odb_filename, static_pickle_filename, cyclic_pickle_filename):
    stress_odb = odbAccess.openOdb(stress_odb_filename)

    instance_name = stress_odb.rootAssembly.instances.keys()[0]
    element_set_names = stress_odb.rootAssembly.instances[instance_name].elementSets.keys()
    element_set_name = None
    for element_set_name in element_set_names:
        if 'ballast_elements' in element_set_name.lower():
            break
    stress_odb.close()
    static_stresses = read_field_from_odb('S', stress_odb_filename, step_name='gravity', set_name=element_set_name,
                                          instance_name=instance_name)
    loading_stresses = read_field_from_odb('S', stress_odb_filename, step_name='loading', set_name=element_set_name,
                                           instance_name=instance_name)

    cyclic_stresses = loading_stresses - static_stresses

    with open(static_pickle_filename, 'wb') as static_pickle:
        pickle.dump({'data': static_stresses, 'instance': instance_name, 'element_set': element_set_name},
                    static_pickle)

    with open(cyclic_pickle_filename, 'wb') as cyclic_pickle:
        pickle.dump({'data': cyclic_stresses, 'instance': instance_name, 'element_set': element_set_name},
                    cyclic_pickle)


if __name__ == '__main__':
    write_stress_pickles(sys.argv[-3], sys.argv[-2], sys.argv[-1])