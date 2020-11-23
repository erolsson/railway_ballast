import sys
import pickle
from abaqus_functions.odb_io_functions import write_field_to_odb


def write_permanent_strains_to_odb(odb_file, pickle_file):
    with open(pickle_file, 'r'):
        data = pickle.load(pickle_file)


if __name__ == '__main__':
    write_permanent_strains_to_odb(sys.argv[-2], sys.argv[-1])