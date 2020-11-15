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
    from abaqusConstants import SINGLE, FIXED, SWEEP, MEDIAL_AXIS, DC3D8, DC3D6, C3D8, C3D6, C3D20, STANDARD, ANALYSIS
    from abaqusConstants import PERCENTAGE, DOMAIN, DEFAULT, INDEX, YZPLANE, XYPLANE, HEX
    from abaqus import backwardCompatibility
    backwardCompatibility.setValues(reportDeprecated=False)
except ImportError:
    print(" ERROR: This script require Abaqus CAE to run")
    raise


class ElasticLayer:
    def __init__(self, name, bottom_width, top_width, height, E, v, density):
        self.name = name
        self.bottom_width = bottom_width
        self.top_width = top_width
        self.height = height
        self.E = E
        self.v = v
        self.density = density


class RailwayEmbankment:
    def __init__(self, layers, length):
        Mdb()
        self.mdb = mdb.models['Model-1']
        self.mdb.setValues(noPartsInputFile=ON)
        self.assembly = self.mdb.rootAssembly
        self.length = length
        self.total_height = 0
        self.layers = layers
        instances = []
        parts = []
        self.sleeper_length = None
        self.sleeper_part = None
        self.sleeper_instances = []
        for layer in layers:
            sketch = self.mdb.ConstrainedSketch(name='sketch_' + layer.name, sheetSize=800.0)
            p1 = (0, self.total_height)
            p2 = (layer.bottom_width, self.total_height)
            p3 = (layer.top_width, self.total_height + layer.height)
            p4 = (0, self.total_height + layer.height)
            self.total_height += layer.height
            sketch.Line(point1=p1, point2=p2)
            sketch.Line(point1=p2, point2=p3)
            sketch.Line(point1=p3, point2=p4)
            sketch.Line(point1=p4, point2=p1)
            part = self.mdb.Part(name='part_' + layer.name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
            part.BaseSolidExtrude(sketch=sketch, depth=length)
            parts.append(part)
            instances.append(self.assembly.Instance(name='instance_' + layer.name, part=part, dependent=ON))
        self.part = self.assembly.PartFromBooleanMerge(name='embankment_part',
                                                       instances=instances,
                                                       keepIntersections=True)
        self.instance = self.assembly.Instance(name='embankment_instance', part=self.part, dependent=ON)
        for instance in instances:
            del self.assembly.instances[instance.name]
        for part in parts:
            del self.mdb.parts[part.name]

        for i in range(len(layers) - 1):
            if layers[i].top_width != layers[i+1].bottom_width:
                datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE,
                                                                            offset=layers[i+1].bottom_width)
                self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                                    cells=self.part.cells)

    def create_sleepers(self, sleeper_cc_distance, sleeper_width, sleeper_height, sleeper_length):
        # Sketching a sleeper
        sketch = self.mdb.ConstrainedSketch(name='sketch_sleeper', sheetSize=800.0)
        p1 = (0, self.total_height)
        p2 = (sleeper_length, self.total_height)
        p3 = (sleeper_length, self.total_height + sleeper_height)
        p4 = (0, self.total_height + sleeper_height)

        sketch.Line(point1=p1, point2=p2)
        sketch.Line(point1=p2, point2=p3)
        sketch.Line(point1=p3, point2=p4)
        sketch.Line(point1=p4, point2=p1)
        self.sleeper_part = part = self.mdb.Part(name='part_sleeper', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.sleeper_part.BaseSolidExtrude(sketch=sketch, depth=sleeper_width)

        # meshing the sleeper
        self.sleeper_part.seedPart(size=0.1)
        elem_type1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
        self.sleeper_part.setElementType(regions=(part.cells,), elemTypes=(elem_type1,))
        self.sleeper_part.generateMesh()

        start_point = sleeper_cc_distance/2 - sleeper_width/2
        sleeper_idx = 0

        datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE,
                                                                    offset=sleeper_length)
        self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                            cells=self.part.cells)
        self.sleeper_length = sleeper_length

        while start_point + sleeper_width < self.length:
            datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                        offset=start_point)
            self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                                cells=self.part.cells)
            datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                        offset=start_point + sleeper_width)
            self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                                cells=self.part.cells)

            instance = self.assembly.Instance(name='sleeper' + str(sleeper_idx), part=part, dependent=ON)
            instance.translate((0, 0, start_point))
            self.sleeper_instances.append(instance)
            sleeper_idx += 1

            # Creating a tie between the sleeper and the ballast
            ballast_instance = self.assembly.instances['embankment_instance']
            ballast_face = ballast_instance.faces.findAt(((1e-3, self.total_height, start_point + 1e-3),))
            ballast_surface = self.assembly.Surface(side1Faces=ballast_face,
                                                    name='ballast_sleeper_surface_' + str(sleeper_idx))
            sleeper_face = instance.faces.findAt(((1e-3, self.total_height, start_point + 1e-3),))
            sleeper_surface = self.assembly.Surface(side1Faces=sleeper_face,
                                                    name='sleeper_surface_' + str(sleeper_idx))
            self.mdb.Tie(name='tie_sleeper_' + str(sleeper_idx), slave=sleeper_surface, master=ballast_surface)
            start_point += sleeper_cc_distance

    def mesh(self):
        def get_edge_direction(edge):
            vertex_idx = edge.getVertices()
            v1 = self.part.vertices[vertex_idx[0]]
            v2 = self.part.vertices[vertex_idx[1]]
            p1 = v1.pointOn[0]
            p2 = v2.pointOn[0]

            dx = p1[0] - p2[0]
            dy = p1[1] - p2[1]
            dz = p1[2] - p2[2]
            return [dx, dy,  dz]

        element_size_ballast_height = 0.1
        element_size_ballast_width = 0.1
        element_size_ballast_length = 0.1
        element_size_subgrade_height = 0.5
        ballast_start_height = 0
        first_ballast_layer = None
        for layer in self.layers:
            first_ballast_layer = layer
            if layer.name.startswith('ballast'):
                break
            ballast_start_height += layer.height

        ballast_top_height = sum([layer.height for layer in self.layers])

        ballast_edges = self.part.edges.getByBoundingBox(xMin=-1e-3,
                                                         xMax=first_ballast_layer.bottom_width + 1e-3,
                                                         yMin=ballast_start_height - 1e-3,
                                                         yMax=ballast_top_height + 1e-3,
                                                         zMin=-1e-3,
                                                         zMax=self.length + 1e-3)
        length_ballast_edges = []
        width_ballast_edges = []
        height_ballast_edges = []
        for ballast_edge in ballast_edges:
            n = get_edge_direction(ballast_edge)
            if n[0] == 0 and n[1] == 0:
                length_ballast_edges.append(ballast_edge)
            if n[1] == 0 and n[2] == 0 and ballast_edge.pointOn[0][0] < self.sleeper_length:
                width_ballast_edges.append(ballast_edge)
            if n[0] == 0 and n[2] == 0 and ballast_edge.pointOn[0][0] < self.sleeper_length:
                height_ballast_edges.append(ballast_edge)

        self.part.seedEdgeBySize(edges=height_ballast_edges, size=element_size_ballast_height, constraint=FIXED)
        self.part.seedEdgeBySize(edges=length_ballast_edges, size=element_size_ballast_length, constraint=FIXED)
        self.part.seedEdgeBySize(edges=width_ballast_edges, size=element_size_ballast_width, constraint=FIXED)
        subgrade_edges = self.part.edges.getByBoundingBox(xMin=-1e-3,
                                                          xMax=self.layers[0].bottom_width + 1e-3,
                                                          yMin=-1e-3, yMax=ballast_start_height + 1e-3,
                                                          zMin=-1e-3, zMax=self.length + 1e-3)
        subgrade_width_edges = []
        subgrade_height_edges = []
        for subgrade_edge in subgrade_edges:
            n = get_edge_direction(subgrade_edge)
            if n[1] == 0 and n[2] == 0 and subgrade_edge.pointOn[0][0] > first_ballast_layer.bottom_width:
                subgrade_width_edges.append(subgrade_edge)
            if n[0] == 0 and n[2] == 0:
                subgrade_height_edges.append(subgrade_edge)

        self.part.seedEdgeBySize(edges=subgrade_height_edges, size=element_size_subgrade_height, constraint=FIXED)
        self.part.seedEdgeBySize(edges=subgrade_width_edges, size=element_size_subgrade_height, constraint=FIXED)

        # Meshing the ballast first
        ballast_cells = self.part.cells.getByBoundingBox(xMax=first_ballast_layer.bottom_width + 1e-3,
                                                         yMin=ballast_start_height - 1e-3)
        self.part.generateMesh(regions=ballast_cells)
        subgrade_cells = self.part.cells.getByBoundingBox(yMax=ballast_start_height + 1e-3)
        self.part.generateMesh(regions=subgrade_cells)
        # elem_type1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
        # self.part.setElementType(regions=(self.part.cells,), elemTypes=(elem_type1,))

    def apply_boundary_conditions(self):
        # Setting displacement boundary conditions
        bottom_instance = self.assembly.instances['embankment_instance']
        bottom_face = bottom_instance.faces.findAt((1e-3, 0, 1e-3)).getFacesByFaceAngle(0.)
        bottom_nodes = self.assembly.Set(name='bottom_nodes', faces=bottom_face)

        x_sym_faces = self.instance.faces.findAt((0, 1e-3, 1e-3)).getFacesByFaceAngle(0.)
        for sleeper in self.sleeper_instances:
            x_sym_faces += sleeper.faces.getByBoundingBox(xMax=1e-3)
        x_sym_nodes = self.assembly.Set(name='x_sym_nodes', faces=x_sym_faces)

        z_sym_faces = self.instance.faces.findAt((1e-3, 1e-3, 0)).getFacesByFaceAngle(0.)
        z_sym_nodes = self.assembly.Set(name='z_sym_nodes', faces=z_sym_faces)

        z1_faces = self.instance.faces.findAt((1e-3, 1e-3, self.length)).getFacesByFaceAngle(0.)
        z1_nodes = self.assembly.Set(name='z1_nodes', faces=z1_faces)

        x1_faces = self.instance.faces.findAt((self.layers[0].bottom_width, 1e-3, 1e-3)).getFacesByFaceAngle(0.)
        x1_nodes = self.assembly.Set(name='x1_nodes', faces=x1_faces)

        self.mdb.DisplacementBC(name='y_bc', region=bottom_nodes, createStepName='Initial', u2=0)
        self.mdb.DisplacementBC(name='x_sym_bc', region=x_sym_nodes, createStepName='Initial', u1=0)
        self.mdb.DisplacementBC(name='z_sym_bc', region=z_sym_nodes, createStepName='Initial', u3=0)
        self.mdb.DisplacementBC(name='x1_bc', region=x1_nodes, createStepName='Initial', u1=0)
        self.mdb.DisplacementBC(name='z1_bc', region=z1_nodes, createStepName='Initial', u3=0)

    def create_materials(self):
        layer_start_height = 0
        for layer in self.layers:
            mat = self.mdb.Material('material_' + layer.name)
            mat.Density(table=((layer.density,),))
            mat.Elastic(table=((layer.E, layer.v),))
            if layer.name.startswith('ballast'):
                mc = mat.MohrCoulombPlasticity(table=((45., 0.),), useTensionCutoff=ON)
                mc.MohrCoulombHardening(table=((1e9, 0.),))
                mc.TensionCutOff(table=((1e4, 0),))

            layer_elements = self.part.elements.getByBoundingBox(yMin=layer_start_height - 1e-3,
                                                                 yMax=layer_start_height + layer.height + 1e-3)
            self.mdb.HomogeneousSolidSection(name='section_' + layer.name,  material='material_' + layer.name)
            self.part.SectionAssignment(region=(layer_elements,), sectionName='section_' + layer.name)
            layer_start_height += layer.height
        mat = self.mdb.Material('concrete')
        mat.Density(table=((2570,),))
        mat.Elastic(table=((30e9, 0.2),))
        self.mdb.HomogeneousSolidSection(name='section_sleeper', material='concrete')
        self.sleeper_part.SectionAssignment(region=(self.sleeper_part.elements,), sectionName='section_sleeper')

    def create_load_steps(self):
        gravity_step = self.mdb.StaticStep(name='gravity', previous='Initial', nlgeom=ON, initialInc=1e-2)
        gravity_step.control.setValues(allowPropagation=OFF,
                                       lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
                                       timeIncrementation=(15, 15, 9.0, 16.0, 10.0, 4.0, 12.0,
                                                           20.0, 6.0, 3.0, 10.0, 0.5, 0.5,
                                                           0.75, 0.85, 0.25, 0.25, 1.5, 0.75),
                                       resetDefaultValues=OFF)
        self.mdb.Gravity(name='gravity', createStepName='gravity', comp2=-9.82)

    def run_job(self, cpus=8):
        job = mdb.Job(name='embankment', model=self.mdb, numCpus=cpus, numDomains=cpus)
        job.submit()
        job.waitForCompletion()


def main():
    embkankment_length = 5.5
    layers = [ElasticLayer(name='soil', bottom_width=20, top_width=20.,
                           height=2., E=10e6, v=0.3, density=2000),
              ElasticLayer(name='clay', bottom_width=20, top_width=20.,
                           height=4., E=7.5e6, v=0.3, density=2000),
              ElasticLayer(name='subgrade', bottom_width=20, top_width=20.,
                           height=1., E=20e6, v=0.35, density=1900),
              ElasticLayer(name='ballast_1', bottom_width=9.3, top_width=3.6,
                           height=3.8, E=200e6, v=0.35, density=1600),
              ElasticLayer(name='ballast_2', bottom_width=2.55, top_width=1.8,
                           height=0.5, E=200e6, v=0.35, density=1600)]
    embankment = RailwayEmbankment(layers, embkankment_length)
    embankment.create_sleepers(sleeper_cc_distance=0.65, sleeper_width=0.265, sleeper_length=1., sleeper_height=0.05)
    embankment.mesh()
    embankment.apply_boundary_conditions()
    embankment.create_materials()
    embankment.create_load_steps()
    embankment.assembly.regenerate()
    embankment.run_job()


if __name__ == '__main__':
    main()
