import pickle
import sys

import numpy as np

from visualization import *
import xyPlot
import odbAccess

from abaqusConstants import POINT_LIST, ELEMENT_NODAL, TRUE_DISTANCE, UNDEFORMED, PATH_POINTS, COMPONENT
from abaqusConstants import NODAL, INTEGRATION_POINT, CENTROID

output_positions = {'ELEMENT_NODAL': ELEMENT_NODAL,
                    'NODAL': NODAL,
                    'INTEGRATION_POINT': INTEGRATION_POINT,
                    'CENTROID': CENTROID}


def create_path(points, path_name, session):
    path_points = []
    for point in points:
        path_points.append((point[0], point[1], point[2]))

    path = session.Path(name=path_name, type=POINT_LIST, expression=path_points)
    return path


def get_data_from_path(path, session, variable, component=None, output_position=ELEMENT_NODAL):
    if component is None:
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variable,
                                                                       outputPosition=output_position)
    else:
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(variableLabel=variable,
                                                                       outputPosition=output_position,
                                                                       refinement=[COMPONENT, component])
    xy = xyPlot.XYDataFromPath(name=path.name + '_' + variable, path=path,
                               labelType=TRUE_DISTANCE, shape=UNDEFORMED, pathStyle=PATH_POINTS,
                               includeIntersections=False)
    return np.array(xy)


def main():
    pickle_file_name = sys.argv[-1]
    with open(pickle_file_name, 'r') as parameter_pickle:
        parameters = pickle.load(parameter_pickle)

    odb_filename = str(parameters['odb_filename'])
    path_points_filename = str(parameters['path_points_filename'])
    variable = str(parameters['variable'])
    output_position = output_positions[str(parameters['output_position'])]
    data_filename = str(parameters['data_filename'])
    component = None
    if 'component' in parameters:
        component = str(parameters['component'])

    odb = odbAccess.openOdb(odb_filename)
    session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=309.913116455078,
                     height=230.809509277344)
    session.viewports['Viewport: 1'].makeCurrent()
    session.viewports['Viewport: 1'].maximize()
    o7 = session.odbs[session.odbs.keys()[0]]
    session.viewports['Viewport: 1'].setValues(displayedObject=o7)

    if 'step_name' not in parameters:
        step_name = odb.steps.keys()[-1]
    else:
        step_name = str(parameters['step_name'])

    step_index = odb.steps.keys().index(step_name)
    if 'frame_number' not in parameters:
        frame_number = len(odb.steps[step_name].frames)
    else:
        frame_number = parameters['frame_number']
    session.viewports['Viewport: 1'].odbDisplay.setFrame(step=step_index, frame=frame_number)

    path_points = np.load(path_points_filename)
    path = create_path(path_points, 'path', session)
    data = get_data_from_path(path, session, variable, component, output_position=output_position)
    np.save(data_filename, data)
    odb.close()


if __name__ == '__main__':
    main()
