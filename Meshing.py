from aircraft_fem.examples.aircraft.geom.glueing import GeneralFuse
from aircraft_fem.examples.aircraft.fem.mesh import MeshGenerator
from parapy.mesh.salome import FixedLength, Mesh as PPMesh, Quad, Tri
from parapy.mesh.salome.grid import SubGrid
from OCC.utils.utilities import make_TopoDS_Compound
from parapy.geom.occ.utilities import topods_shape_getter
from torsionbox import TorsionBox
from parapy.core import *

# Add these imports at the top of torsionbox.py if not already present
from parapy.core import Attribute, Part


class Mesh(PPMesh):
    def get_subgrid_on_the_fly(self, shape, *other_shapes, label):
        if other_shapes:
            topods_shape = make_TopoDS_Compound(map(
                topods_shape_getter, [shape, *other_shapes]))
        else:
            topods_shape = shape.TopoDS_Shape
        smesh_mesh = self.SMESH_Mesh
        smesh_submesh = smesh_mesh.GetSubMesh(topods_shape)
        return SubGrid(SMESH_Mesh=smesh_mesh,
                       SMESH_subMesh=smesh_submesh,
                       dimension=shape.TOPODIM,
                       label=label)


class MeshGenerator(Base):
    shape_to_mesh = Input()
    element_length = Input(0.2)

    @Attribute
    def quad_faces(self):
        faces = [f for f in self.shape_to_mesh.faces if len(f.edges) == 4]
        lst = []

        taper_ratio = 0.2
        for f in faces:
            e1 = f.edges[0]
            try:
                e3 = e1.opposite_edge
            except Exception:
                continue

            if abs(e1.length - e3.length) > taper_ratio * e1.length:
                continue

            e2 = f.edges[1]
            e4 = e2.opposite_edge
            if abs(e2.length - e4.length) > taper_ratio * e2.length:
                continue

            lst.append(f)
        return lst

    @Part(in_tree=False)
    def fixed_length(self):
        return FixedLength(shape_to_mesh=self.shape_to_mesh,
                           length=self.element_length)

    @Part(in_tree=False)
    def quad(self):
        return Quad(quantify=len(self.quad_faces),
                    shape=self.quad_faces[child.index])

    @Part
    def tri(self):
        return Tri(shape_to_mesh=self.shape_to_mesh,
                   quad_dominant=False,
                   only_2d=True,
                   min_size=0.1,
                   max_size=0.3)

    @Part
    def mesh(self):
        return Mesh(shape_to_mesh=self.shape_to_mesh,
                    display_mode="shaded",
                    controls=[self.fixed_length, self.quad, self.tri])


class FinalMesh(Base):
    @Attribute
    def torsionbox(self):
        return TorsionBox()

    @Part
    def shape_to_mesh(self):
        return GeneralFuse(
            tools=[
                      self.torsionbox.wing_upper_surface,
                      self.torsionbox.wing_lower_surface,
                      self.torsionbox.front_spar,
                      self.torsionbox.rear_spar
                  ] + [rib.rib_surface for rib in self.torsionbox.ribs],

            transparency=0.5
        )

    @Part
    def mesh_generator(self):
        return MeshGenerator(shape_to_mesh=self.shape_to_mesh)

    @Attribute
    def mesh(self):
        return self.mesh_generator.mesh


if __name__ == '__main__':
    from parapy.gui import display

    obj = FinalMesh(label="Mesh Torsion Box")
    display(obj)
