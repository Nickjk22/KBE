# Import necessary modules and classes
from glueing import GeneralFuse
from parapy.mesh.salome import FixedLength, Mesh as PPMesh, Quad, Tri
from parapy.mesh.salome.grid import SubGrid
from OCC.utils.utilities import make_TopoDS_Compound
from parapy.geom.occ.utilities import topods_shape_getter
from torsionbox import TorsionBox
from parapy.core import *
from parapy.core import Attribute, Part

# Extended Mesh class with a helper to get subgrids for selected shapes
class Mesh(PPMesh):
    def get_subgrid_on_the_fly(self, shape, *other_shapes, label):
        # Create compound shape if multiple shapes are provided
        if other_shapes:
            topods_shape = make_TopoDS_Compound(map(
                topods_shape_getter, [shape, *other_shapes]))
        else:
            topods_shape = shape.TopoDS_Shape
        # Get the mesh and submesh from SMESH API
        smesh_mesh = self.SMESH_Mesh
        smesh_submesh = smesh_mesh.GetSubMesh(topods_shape)
        return SubGrid(SMESH_Mesh=smesh_mesh,
                       SMESH_subMesh=smesh_submesh,
                       dimension=shape.TOPODIM,
                       label=label)

# MeshGenerator class: generates wing mesh from parametric input
class MeshGenerator(Base):
    # Geometry and control inputs
    check_element = Input()

    wing_airfoil_root = Input()
    wing_airfoil_middle = Input()
    wing_airfoil_tip = Input()

    wing_root_chord = Input()
    wing_middle_chord = Input()
    wing_tip_chord = Input()

    wing_thickness_factor_root = Input()
    wing_thickness_factor_middle = Input()
    wing_thickness_factor_tip = Input()

    wing_semi_span_planform1 = Input()
    wing_semi_span = Input()
    wing_sweep_leading_edge_planform1 = Input()
    wing_sweep_leading_edge_planform2 = Input()
    wing_twist = Input()

    # Spar positions
    front_spar_position = Input()
    rear_spar_position = Input()

    # Number of ribs
    rib_number = Input()

    # Discretization settings
    section_number = Input()
    points_number = Input()
    element_length = Input()

    @Attribute
    def torsionbox(self):
        # Generates the 3D wing structure using a torsion box representation
        return TorsionBox(wing_airfoil_root=self.wing_airfoil_root,
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

                          front_spar_position=self.front_spar_position,
                          rear_spar_position=self.rear_spar_position,
                          rib_number=self.rib_number,

                          section_number=self.section_number,
                          points_number=self.points_number,
                          hidden=True)

    @Part
    def shape_to_mesh(self):
        # Merge and fuse all structural components into one compound for meshing
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

    @Attribute
    def quad_faces(self):
        # Select faces that are suitable for quadrilateral meshing
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
        # Fixed edge length control for mesh generation
        return FixedLength(shape_to_mesh=self.shape_to_mesh,
                           length=self.element_length)

    @Part(in_tree=False)
    def quad(self):
        # Quadrilateral mesh generation for selected faces
        return Quad(quantify=len(self.quad_faces),
                    shape=self.quad_faces[child.index])

    @Part
    def tri(self):
        # Triangular mesh for all surfaces (fallback or full domain)
        return Tri(shape_to_mesh=self.shape_to_mesh,
                   quad_dominant=False,
                   only_2d=False,
                   min_size=0.01,
                   max_size=0.3,
                   hidden=True)

    @Part
    def mesh(self):
        # Full mesh object combining fixed length and triangulation
        return Mesh(shape_to_mesh=self.shape_to_mesh,
                    display_mode="shaded",
                    controls=[self.fixed_length, self.tri])

    @Attribute
    def face_hash_map(self):
        # Create a hash map for mesh faces to track identity/changes
        return [(i, hash(face.TopoDS_Shape)) for i, face in enumerate(self.shape_to_mesh.faces)]


# Standalone visualization if this script is run directly
if __name__ == '__main__':
    from parapy.gui import display

    obj = MeshGenerator(label="Mesh Torsion Box",
                        check_element=0,
                        wing_airfoil_root="whitcomb_interpolated.dat",
                        wing_airfoil_middle="whitcomb_interpolated.dat",
                        wing_airfoil_tip="whitcomb_interpolated.dat",

                        wing_root_chord=6,
                        wing_middle_chord=4,
                        wing_tip_chord=1.5,

                        wing_thickness_factor_root=1,
                        wing_thickness_factor_middle=1,
                        wing_thickness_factor_tip=1,

                        wing_semi_span_planform1=5,
                        wing_semi_span=16,
                        wing_sweep_leading_edge_planform1=20,
                        wing_sweep_leading_edge_planform2=20,
                        wing_twist=0,

                        front_spar_position=0.2,
                        rear_spar_position=0.6,
                        rib_number=12,

                        section_number=14,
                        points_number=14,
                        element_length=0.1
                        )

    display(obj)
