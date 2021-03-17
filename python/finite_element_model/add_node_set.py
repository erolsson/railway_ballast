import os

import odbAccess


def add_node_set_to_odb(odb_file_name, node_set_name, x_min=-1e99, x_max=1e99, y_min=-1e99, y_max=1e99,
                        z_min=-1e99, z_max=1e99, instance_name=None):
    odb = odbAccess.openOdb(odb_file_name, readOnly=False)
    if instance_name is None:
        instance_name = odb.rootAssembly.instances.keys()[0]
    nodes = odb.rootAssembly.instances[instance_name].nodes
    set_node_labels = []
    for node in nodes:
        x, y, z = node.coordinates
        if x_min < x < x_max and y_min < y < y_max and z_min < z < z_max:
            set_node_labels.append(node.label)

    odb.rootAssembly.instances[instance_name].NodeSetFromNodeLabels(name=node_set_name, nodeLabels=set_node_labels)
    odb.save()
    odb.close()


if __name__ == '__main__':
    odb_directory = os.path.expanduser('~/railway_ballast/odbs/')
    add_node_set_to_odb(odb_directory + 'embankment_sleepers_low_30_0t.odb', 'ballast_bottom_nodes',
                        y_min=7-1e-3, y_max=7+1e-3)
    add_node_set_to_odb(odb_directory + 'embankment_sleepers_high_30_0t.odb',
                        'ballast_bottom_nodes', y_min=7-1e-3, y_max=7+1e-3)
