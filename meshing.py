from aircraft_fem.examples.aircraft.geom.glueing import GeneralFuse
from parapy.mesh.salome import FixedLength, Mesh as PPMesh, Quad, Tri
from parapy.mesh.salome.grid import SubGrid
from OCC.utils.utilities import make_TopoDS_Compound
from parapy.geom.occ.utilities import topods_shape_getter
from torsionbox import TorsionBox
from parapy.core import *
from parapy.geom.occ.brep import BRep
from sections import Section
from parapy.mesh import EdgeGroup, FaceGroup

import numpy as np

# Add these imports at the top of torsionbox.py if not already present
from parapy.core import Attribute, Part

FACE = "face_group"
CONSTRAINED_EDGE1 = "constrained_edge1_group"
CONSTRAINED_EDGE2 = "constrained_edge2_group"
LOADED_EDGE = "loaded_edge_group"

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
    element_length = Input(0.5)
    wing_span = Input(30)
    sections = Input(14)

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
                   only_2d=False,
                   min_size=0.1,
                   max_size=0.3,
                   hidden=True)

    @Part
    def mesh(self):
        return Mesh(shape_to_mesh=self.shape_to_mesh,
                    display_mode="shaded",
                    controls=[self.fixed_length, self.tri])


class FinalMesh(Base):
    check_element = Input(0)

    section_number = Input(14)

    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

    wing_root_chord = Input(6)
    wing_middle_chord = Input(4)
    wing_tip_chord = Input(1.5)

    wing_thickness_factor_root = Input(1)
    wing_thickness_factor_middle = Input(1)
    wing_thickness_factor_tip = Input(1)

    wing_semi_span_planform1 = Input(5)
    wing_semi_span = Input(16)
    wing_sweep_leading_edge_planform1 = Input(20)
    wing_sweep_leading_edge_planform2 = Input(20)
    wing_twist = Input(0)

    # Spars
    front_spar_thickness = Input(1)
    front_spar_position = Input(0.2)
    rear_spar_thickness = Input(1)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_thickness = Input(0.2)
    rib_number = Input(12)

    @Part
    def sections(self):
        return Section(wing_airfoil_root=self.wing_airfoil_root,
                       wing_airfoil_middle=self.wing_airfoil_middle,
                       wing_airfoil_tip=self.wing_airfoil_tip,

                       wing_root_chord=self.wing_root_chord,
                       wing_middle_chord=self.wing_middle_chord,
                       wing_tip_chord=self.wing_tip_chord,

                       wing_thickness_factor_root=self.wing_thickness_factor_root,
                       wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                       wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                       wing_semi_span_planform1=self.wing_semi_span_planform1,
                       wing_semi_span=self.wing_semi_span,
                       wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                       wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,
                       wing_twist=self.wing_twist,

                       section_number=self.section_number,
                       hidden=True
                       )

    @Attribute
    def torsionbox(self):
        return TorsionBox()

    @Part
    def shape_to_mesh(self):
        return GeneralFuse(
            tools=[self.torsionbox.wing_upper_surface,
                   self.torsionbox.wing_lower_surface,
                   self.torsionbox.front_spar,
                   self.torsionbox.rear_spar
                   ] + [rib.rib_surface for rib in self.torsionbox.ribs],
            transparency=0.9,
            mesh_deflection=0.001,
            fuzzy_value=0.001
        )

    @Part
    def mesh_generator(self):
        return MeshGenerator(shape_to_mesh=self.shape_to_mesh)

    @Attribute
    def mesh(self):
        return self.mesh_generator.mesh

    @Part
    def highlighted_face(self):
        return BRep(TopoDS_Shape=self.shape_to_mesh.faces[self.check_element].TopoDS_Shape,
                    color="red",
                    transparency=0)

    @Attribute
    def number_of_faces(self):
        return len(self.shape_to_mesh.faces)

    @Attribute
    def test(self):
        return self.shape_to_mesh.edges

    @Attribute
    def face_hash_map(self):
        return [(i, hash(face.TopoDS_Shape)) for i, face in enumerate(self.shape_to_mesh.faces)]


if __name__ == '__main__':
    from parapy.gui import display

    obj = FinalMesh(label="Mesh Torsion Box")

    display(obj)
