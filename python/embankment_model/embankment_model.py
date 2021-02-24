from __future__ import print_function, division

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
    from abaqusConstants import PERCENTAGE, DOMAIN, DEFAULT, INDEX, YZPLANE, XYPLANE, HEX, TOTAL_FORCE, NUMBER
    from abaqus import backwardCompatibility
    backwardCompatibility.setValues(reportDeprecated=False)
except ImportError:
    print(" ERROR: This script require Abaqus CAE to run")
    raise

from simulations import embankment_21_22_5t as simulation_to_run


class RailwayEmbankment:
    def __init__(self, simulation):
        Mdb()
        self.mdb = mdb.models['Model-1']
        self.mdb.setValues(noPartsInputFile=ON)
        self.assembly = self.mdb.rootAssembly
        self.length = simulation.embankment_length
        self.total_height = 0
        self.layers = simulation.layers
        instances = []
        parts = []
        self.sleeper_length = None
        self.sleeper_part = None
        self.center_sleeper_part = None
        self.sleeper_instances = []
        self.rail_part = None
        self.rail_instance = None
        self.track_gauge = simulation.track_gauge
        self.sleeper_height = 0.175

        self.rail_area_inertia = 32520000/1e12
        self.rail_height = 0.172
        self.rail_width = 12*self.rail_area_inertia/self.rail_height**3
        self.rail_density = 60.4/self.rail_width/self.rail_height
        self.job_name = simulation.job_name

        for layer in self.layers:
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
            part.BaseSolidExtrude(sketch=sketch, depth=self.length)
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

        for i in range(len(self.layers) - 1):
            if self.layers[i].top_width != self.layers[i+1].bottom_width:
                datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE,
                                                                            offset=self.layers[i+1].bottom_width)
                self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                                    cells=self.part.cells)

    def create_sleepers(self, sleeper_cc_distance, sleeper_width, sleeper_length, center_sleeper=True):
        # Sketching a sleeper
        sketch = self.mdb.ConstrainedSketch(name='sketch_sleeper', sheetSize=800.0)
        p1 = (0, self.total_height)
        p2 = (sleeper_length, self.total_height)
        p3 = (sleeper_length, self.total_height + self.sleeper_height)
        p4 = (0, self.total_height + self.sleeper_height)

        sketch.Line(point1=p1, point2=p2)
        sketch.Line(point1=p2, point2=p3)
        sketch.Line(point1=p3, point2=p4)
        sketch.Line(point1=p4, point2=p1)
        sleeper_names = ['sleeper']
        sleeper_widths = [sleeper_width]
        if center_sleeper:
            sleeper_names.append('center_sleeper')
            sleeper_widths.append(sleeper_width/2)
        for name, width in zip(sleeper_names, sleeper_widths):
            part = self.mdb.Part(name='part_' + name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
            part.BaseSolidExtrude(sketch=sketch, depth=width)

            x = self.track_gauge/2 - self.rail_width/2
            datum_plane_vertical = part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=x)
            part.PartitionCellByDatumPlane(datumPlane=part.datum[datum_plane_vertical.id], cells=part.cells)

            x = self.track_gauge/2 + self.rail_width/2
            datum_plane_vertical = part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE, offset=x)
            part.PartitionCellByDatumPlane(datumPlane=part.datum[datum_plane_vertical.id], cells=part.cells)

            # meshing the sleeper
            part.seedPart(size=0.025)
            elem_type1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
            part.setElementType(regions=(part.cells,), elemTypes=(elem_type1,))
            part.generateMesh()
            if name == 'sleeper':
                self.sleeper_part = part
            else:
                self.center_sleeper_part = part

        datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=YZPLANE,
                                                                    offset=sleeper_length)
        self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                            cells=self.part.cells)
        self.sleeper_length = sleeper_length
        if center_sleeper:
            start_point = sleeper_cc_distance - sleeper_width/2
            datum_plane_vertical = self.part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                        offset=sleeper_width/2)
            self.part.PartitionCellByDatumPlane(datumPlane=self.part.datum[datum_plane_vertical.id],
                                                cells=self.part.cells)
            datum_plane_vertical = self.rail_part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                             offset=sleeper_width/2)
            self.rail_part.PartitionCellByDatumPlane(datumPlane=self.rail_part.datum[datum_plane_vertical.id],
                                                     cells=self.rail_part.cells)
            instance = self.assembly.Instance(name='sleeper_0', part=self.center_sleeper_part, dependent=ON)
            self.sleeper_instances.append(instance)
            ballast_instance = self.assembly.instances['embankment_instance']
            ballast_face = ballast_instance.faces.findAt(((1e-3, self.total_height, 1e-3),))
            ballast_surface = self.assembly.Surface(side1Faces=ballast_face,
                                                    name='ballast_sleeper_surface_0')
            sleeper_face = instance.faces.findAt((1e-3, self.total_height, 1e-3)).getFacesByFaceAngle(0.)
            sleeper_surface = self.assembly.Surface(side1Faces=sleeper_face,
                                                    name='lower_sleeper_surface_0')
            self.mdb.Tie(name='tie_sleeper_0', slave=sleeper_surface, master=ballast_surface)

            sleeper_face = instance.faces.findAt(((self.track_gauge/2,
                                                   self.total_height + self.sleeper_height,
                                                   1e-3),))
            sleeper_surface = self.assembly.Surface(side1Faces=sleeper_face,
                                                    name='upper_sleeper_surface_0')
            rail_face = self.rail_instance.faces.findAt(((self.track_gauge/2,
                                                          self.total_height + self.sleeper_height,
                                                          1e-3),))
            rail_surface = self.assembly.Surface(side1Faces=rail_face,
                                                 name='rail_lower_surface_0')
            self.mdb.Tie(name='tie_rail_0', slave=sleeper_surface, master=rail_surface)

        else:
            start_point = sleeper_cc_distance/2 - sleeper_width/2
        sleeper_idx = 1

        while start_point + sleeper_width < self.length:
            for part in [self.part, self.rail_part]:
                datum_plane_vertical = part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                       offset=start_point)
                part.PartitionCellByDatumPlane(datumPlane=part.datum[datum_plane_vertical.id], cells=part.cells)
                datum_plane_vertical = part.DatumPlaneByPrincipalPlane(principalPlane=XYPLANE,
                                                                       offset=start_point + sleeper_width)
                part.PartitionCellByDatumPlane(datumPlane=part.datum[datum_plane_vertical.id], cells=part.cells)

            instance = self.assembly.Instance(name='sleeper' + str(sleeper_idx), part=self.sleeper_part, dependent=ON)
            instance.translate((0, 0, start_point))
            self.sleeper_instances.append(instance)

            # Creating a tie between the sleeper and the ballast and between the rail and the sleeper
            ballast_instance = self.assembly.instances['embankment_instance']
            ballast_face = ballast_instance.faces.findAt(((1e-3, self.total_height, start_point + 1e-3),))
            ballast_surface = self.assembly.Surface(side1Faces=ballast_face,
                                                    name='ballast_sleeper_surface_' + str(sleeper_idx))
            sleeper_face = instance.faces.findAt((1e-3, self.total_height, start_point + 1e-3)).getFacesByFaceAngle(0.)
            sleeper_surface = self.assembly.Surface(side1Faces=sleeper_face,
                                                    name='sleeper_surface_' + str(sleeper_idx))
            self.mdb.Tie(name='tie_sleeper_' + str(sleeper_idx), slave=ballast_surface, master=sleeper_surface)

            rail_face = self.rail_instance.faces.findAt(((self.track_gauge/2,
                                                         self.total_height + self.sleeper_height,
                                                         start_point + 1e-3),))
            sleeper_face = instance.faces.findAt(((self.track_gauge/2,
                                                   self.total_height + self.sleeper_height,
                                                   start_point + 1e-3),))
            sleeper_surface = self.assembly.Surface(side1Faces=(sleeper_face,),
                                                    name='upper_sleeper_surface_' + str(sleeper_idx))
            rail_surface = self.assembly.Surface(side1Faces=rail_face,
                                                 name='rail_surface_' + str(sleeper_idx))
            self.mdb.Tie(name='tie_rail_' + str(sleeper_idx), slave=sleeper_surface, master=rail_surface)
            start_point += sleeper_cc_distance

            sleeper_idx += 1

    def create_rail(self):
        sketch = self.mdb.ConstrainedSketch(name='sketch_sleeper', sheetSize=800.0)
        p1 = (self.track_gauge/2 - self.rail_width/2, self.total_height + self.sleeper_height)
        p2 = (self.track_gauge/2 + self.rail_width/2, self.total_height + self.sleeper_height)
        p3 = (self.track_gauge/2 + self.rail_width/2, self.total_height + self.sleeper_height + self.rail_height)
        p4 = (self.track_gauge/2 - self.rail_width/2, self.total_height + self.sleeper_height + self.rail_height)

        sketch.Line(point1=p1, point2=p2)
        sketch.Line(point1=p2, point2=p3)
        sketch.Line(point1=p3, point2=p4)
        sketch.Line(point1=p4, point2=p1)

        self.rail_part = self.mdb.Part(name='part_rail', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        self.rail_part.BaseSolidExtrude(sketch=sketch, depth=self.length)
        self.rail_instance = self.assembly.Instance(name='rail_instance', part=self.rail_part, dependent=ON)

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
        element_size_subgrade_width = 0.5
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
        width_ballast_edges_size_seeds = []
        height_ballast_edges = []
        seed_layer_y = self.total_height - self.layers[-1].height
        for ballast_edge in ballast_edges:
            n = get_edge_direction(ballast_edge)
            if n[0] == 0 and n[1] == 0:
                length_ballast_edges.append(ballast_edge)

            if n[1] == 0 and n[2] == 0 and ballast_edge.pointOn[0][1] == seed_layer_y:
                width_ballast_edges_size_seeds.append(ballast_edge)
            if n[0] == 0 and n[2] == 0 and ballast_edge.pointOn[0][0] < self.sleeper_length:
                height_ballast_edges.append(ballast_edge)

        self.part.seedEdgeBySize(edges=height_ballast_edges, size=element_size_ballast_height, constraint=FIXED)
        self.part.seedEdgeBySize(edges=length_ballast_edges, size=element_size_ballast_length, constraint=FIXED)
        self.part.seedEdgeBySize(edges=width_ballast_edges_size_seeds, size=element_size_ballast_width,
                                 constraint=FIXED)

        vertical_edges = {}
        y0 = self.total_height
        for layer in self.layers[::-1]:
            y0 -= layer.height
            y_point = (y0 + y0 + layer.height)/2
            vertical_edges[y_point] = self.part.edges.findAt((0., y_point, 0))

        width_edges = {}
        width_edges_list = self.part.edges.getByBoundingBox(yMin=seed_layer_y - 1e-3, yMax=seed_layer_y + 1e-3,
                                                            zMin=-1e-3, zMax=1e-3)
        for width_edge in width_edges_list:
            vertex_idx = width_edge.getVertices()
            p1 = self.part.vertices[vertex_idx[0]].pointOn[0]
            p2 = self.part.vertices[vertex_idx[1]].pointOn[0]
            x0 = min(p1[0], p2[0])
            width_edges[x0] = width_edge

        for ballast_edge in ballast_edges:
            n = get_edge_direction(ballast_edge)
            if n[1] != 0:
                vertex_idx = ballast_edge.getVertices()
                p1 = self.part.vertices[vertex_idx[0]].pointOn[0]
                p2 = self.part.vertices[vertex_idx[1]].pointOn[0]
                y = (p1[1] + p2[1])/2
                n_seeds = self.part.getEdgeSeeds(edge=vertical_edges[y], attribute=NUMBER)
                self.part.seedEdgeByNumber(edges=[ballast_edge], number=n_seeds)

        for edge in self.part.edges.getByBoundingBox(xMax=first_ballast_layer.bottom_width + 1e-3):
            n = get_edge_direction(edge)
            if n[1] == 0 and n[2] == 0:
                vertex_idx = edge.getVertices()
                p1 = self.part.vertices[vertex_idx[0]].pointOn[0]
                p2 = self.part.vertices[vertex_idx[1]].pointOn[0]
                x0 = min(p1[0], p2[0])
                n_seeds = self.part.getEdgeSeeds(edge=width_edges[x0], attribute=NUMBER)
                self.part.seedEdgeByNumber(edges=[edge], number=n_seeds)

        subgrade_edges = self.part.edges.getByBoundingBox(xMin=-self.layers[0].bottom_width + 1e-3,
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
        self.part.seedEdgeBySize(edges=subgrade_width_edges, size=element_size_subgrade_width, constraint=FIXED)

        # Meshing the ballast first
        ballast_cells = self.part.cells.getByBoundingBox(xMax=first_ballast_layer.bottom_width + 1e-3,
                                                         yMin=ballast_start_height - 1e-3)
        self.part.generateMesh(regions=ballast_cells)
        subgrade_cells = self.part.cells.getByBoundingBox(yMax=ballast_start_height + 1e-3)
        self.part.generateMesh(regions=subgrade_cells)

        if self.rail_part:
            self.rail_part.seedPart(size=0.05)
            elem_type1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
            self.rail_part.setElementType(regions=(self.rail_part.cells,), elemTypes=(elem_type1,))
            self.rail_part.generateMesh()
        elem_type1 = mesh.ElemType(elemCode=C3D20, elemLibrary=STANDARD)
        self.part.setElementType(regions=(self.part.cells,), elemTypes=(elem_type1,))

        ballast_elements = self.part.elements.getByBoundingBox(yMin=ballast_start_height - 1e-3,
                                                               yMax=self.total_height + 1e-3)
        self.part.Set(name='ballast_elements', elements=ballast_elements)
        ballast_nodes = self.part.nodes.getByBoundingBox(xMax=first_ballast_layer.bottom_width + 1e-3,
                                                         yMin=ballast_start_height - 1e-3,
                                                         yMax=self.total_height + 1e-3)
        self.part.Set(name='ballast_nodes', elements=ballast_nodes)

        top_node = self.part.nodes.getByBoundingBox(xMax=1e-3, yMin=self.total_height - 1e-3, zMax=1e-3)
        self.part.Set(name='top_node', nodes=top_node)

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
        z1_faces = self.instance.faces.findAt((1e-3, 1e-3, self.length)).getFacesByFaceAngle(0.)
        if self.center_sleeper_part:
            z_sym_faces += self.sleeper_instances[0].faces.getByBoundingBox(zMax=1e-3)
        if self.rail_part:
            z_sym_faces += self.rail_instance.faces.getByBoundingBox(zMax=1e-3)
            z1_faces += self.rail_instance.faces.getByBoundingBox(zMin=self.length - 1e-3)
        z_sym_nodes = self.assembly.Set(name='z_sym_nodes', faces=z_sym_faces)

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

        if self.sleeper_instances:
            mat = self.mdb.Material('concrete')
            mat.Density(table=((2570,),))
            mat.Elastic(table=((30e9, 0.2),))
            self.mdb.HomogeneousSolidSection(name='section_sleeper', material='concrete')
            self.sleeper_part.SectionAssignment(region=(self.sleeper_part.elements,), sectionName='section_sleeper')
            if self.center_sleeper_part:
                self.center_sleeper_part.SectionAssignment(region=(self.center_sleeper_part.elements,),
                                                           sectionName='section_sleeper')

        if self.rail_instance:
            mat = self.mdb.Material('rail')
            mat.Density(table=((self.rail_density,),))
            mat.Elastic(table=((200e9, 0.3),))
            self.mdb.HomogeneousSolidSection(name='section_rail', material='rail')
            self.rail_part.SectionAssignment(region=(self.rail_part.elements,), sectionName='section_rail')

    def create_load_steps(self, axes_load):
        gravity_step = self.mdb.StaticStep(name='gravity', previous='Initial', nlgeom=ON, initialInc=1e-2,
                                           maxNumInc=10000000)
        gravity_step.control.setValues(allowPropagation=OFF,
                                       lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
                                       timeIncrementation=(15, 15, 9.0, 16.0, 10.0, 4.0, 12.0,
                                                           20.0, 6.0, 3.0, 10.0, 0.5, 0.5,
                                                           0.75, 0.85, 0.25, 0.25, 1.5, 0.75),
                                       resetDefaultValues=OFF)
        self.mdb.Gravity(name='gravity', createStepName='gravity', comp2=-9.82)

        loading_step = self.mdb.StaticStep(name='loading', previous='gravity', nlgeom=ON, initialInc=1e-2,
                                           maxNumInc=10000000)
        loading_step.control.setValues(allowPropagation=OFF,
                                       lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
                                       timeIncrementation=(15, 15, 9.0, 16.0, 10.0, 4.0, 12.0,
                                                           20.0, 6.0, 3.0, 10.0, 0.5, 0.5,
                                                           0.75, 0.85, 0.25, 0.25, 1.5, 0.75),
                                       resetDefaultValues=OFF)
        load_faces = self.rail_instance.faces.findAt(((self.track_gauge/2,
                                                       self.total_height + self.sleeper_height + self.rail_height,
                                                       1e-3),))
        axes_force = axes_load*1000*9.82/4
        load_region = self.assembly.Surface(side1Faces=load_faces, name='load_face')
        self.mdb.Pressure(name='train_wheel', createStepName='loading', region=load_region, magnitude=axes_force,
                          distributionType=TOTAL_FORCE)

    def run_job(self, cpus=12):
        job = mdb.Job(name=self.job_name, model=self.mdb, numCpus=cpus, numDomains=cpus)
        job.submit()
        job.waitForCompletion()


def main():
    embankment = RailwayEmbankment(simulation_to_run)
    embankment.create_rail()
    embankment.create_sleepers(sleeper_cc_distance=0.65, sleeper_width=0.265, sleeper_length=1.)
    embankment.mesh()
    embankment.apply_boundary_conditions()
    embankment.create_materials()
    embankment.create_load_steps(axes_load=simulation_to_run.axes_load)
    embankment.assembly.regenerate()
    embankment.run_job()


if __name__ == '__main__':
    main()
