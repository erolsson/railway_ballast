from __future__ import print_function, division

import sys
import os
from math import sqrt, sin, cos, pi
try:
    import mesh
    from abaqus import session, Mdb, mdb
    import assembly
    import section
    import step
    import mesh
    import interaction
    from abaqusConstants import COORDINATE, STANDALONE, ON, DEFORMABLE_BODY, AXISYM, OFF, THREE_D, DELETE, GEOMETRY
    from abaqusConstants import SINGLE, FIXED, SWEEP, MEDIAL_AXIS, DC3D8, DC3D6, C3D8, C3D6, STANDARD, ANALYSIS
    from abaqusConstants import PERCENTAGE, DOMAIN, DEFAULT, INDEX, YZPLANE, XYPLANE
    from abaqus import backwardCompatibility
    backwardCompatibility.setValues(reportDeprecated=False)
except ImportError:
    print(" ERROR: This script require Abaqus CAE to run")
    raise


class ElasticLayer:
    def __init__(self, name, bottom_width, top_with, height, E, v):
        self.name = name
        self.bottom_width = bottom_width
        self.top_width = top_with
        self.height = height
        self.E = E
        self.v = v


embkankment_length = 40
layers = {}

