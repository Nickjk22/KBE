# --- ParaPy and external library imports ---
from parapy.core import *
from thickness import DEFAULT_SKIN_THICKNESS, THICKNESSES, THICKNESS_TO_RGB
from parapy.core.widgets import Dropdown
from material import STEEL
from parapy.gui.events import EVT_SELECTION_CHANGING
from parapy.geom import ModifiedShape
from parapy.core.decorators import on_event
from torsionbox import TorsionBox
from meshing import MeshGenerator
import numpy as np


# --- Skin class: represents a shell with customizable thickness and color ---
class Skin(Base):
    wing_shell = Input()  # Full wing shell geometry
    faces_to_keep = Input()  # Subset of faces to preserve
    thickness = Input(DEFAULT_SKIN_THICKNESS, widget=Dropdown(THICKNESSES))  # Thickness with dropdown widget
    material = Input(STEEL)  # Default material is steel

    @on_event(EVT_SELECTION_CHANGING)
    def handle_selection_clicks(self, evt):
        # Captures selection events on this object and its subshapes
        evt.owner = self

    @Part
    def shell(self):
        # Creates a shell based on selected faces with color mapped to thickness
        return ModifiedShape(self.wing_shell,
                             keep=self.faces_to_keep,
                             color=THICKNESS_TO_RGB[self.thickness])


# --- CodeAster_primitives class: defines structure and load locations for FEM ---
class CodeAster_primitives(Base):
    # Airfoil inputs (default is Whitcomb)
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

    # Chord lengths
    wing_root_chord = Input(6)
    wing_middle_chord = Input(4)
    wing_tip_chord = Input(1.5)

    # Thickness scaling
    wing_thickness_factor_root = Input(1)
    wing_thickness_factor_middle = Input(1)
    wing_thickness_factor_tip = Input(1)

    # Geometry layout
    wing_semi_span_planform1 = Input(5)
    wing_semi_span = Input(16)
    wing_sweep_leading_edge_planform1 = Input(20)
    wing_sweep_leading_edge_planform2 = Input(20)

    # Spar locations
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    # Ribs and resolution
    rib_number = Input(12)
    section_number = Input(14)
    points_number = Input(14)
    element_length = Input(0.1)

    @Input
    def torsionbox(self):
        # Constructs the internal structure of the wing
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

    @Input
    def finalmesh(self):
        # Generates mesh using the TorsionBox geometry
        return MeshGenerator(check_element=0,
                             wing_airfoil_root=self.wing_airfoil_root,
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
                             element_length=self.element_length)

    @Part
    def skin(self):
        # Defines the structural skin surface for Code_Aster analysis
        return Skin(faces_to_keep=self.torsionbox.wing_upper_surface,
                    wing_shell=self.finalmesh.shape_to_mesh)

    @Attribute
    def spanwise_points_list(self):
        # Normalized list of points along the span (0 to 1)
        return np.linspace(0, 1, self.points_number)

    @Attribute
    def points_list(self):
        # Actual point coordinates from the torsionbox
        return [pt.point for pt in self.torsionbox.points]

    @Attribute
    def load_primitives(self):
        # Finds nearest mesh nodes to specific force application points
        lst = []
        mesh = self.finalmesh.mesh
        for pt in self.points_list:
            tolerance = self.finalmesh.element_length * 2
            node = mesh.grid.find_node_at(pt, tolerance=tolerance)
            if node is None:
                msg = "No mesh node found near strip force cp point. Increase" \
                      " the tolerance of find_node_at."
                raise RuntimeError(msg)
            lst.append(node)
        return lst

    @Attribute
    def foo(self):
        # Reconstructs subgrids from shape history for the skin faces
        lst = []
        history = self.finalmesh.shape_to_mesh.history
        mesh = self.finalmesh.mesh
        get_subgrid = mesh.get_subgrid_on_the_fly

        for face in self.skin.shell.faces:
            shape = history(face)  # gets underlying shape(s)
            label = 'group_ma_' + str(id(face))
            subgrid = get_subgrid(shape[0], *shape[1:], label=label)
            lst.append(subgrid)

        return lst

    @Attribute
    def structural_elements(self):
        # All major load-bearing structural parts of the torsionbox
        return (
                self.torsionbox.ribs +
                [self.torsionbox.front_spar] +
                [self.torsionbox.rear_spar] +
                [self.torsionbox.wing_upper_surface] +
                [self.torsionbox.wing_lower_surface]
        )

    @property
    def structural_element_primitives(self):
        # Flatten all structural faces (if applicable) into a single list
        elements = []
        for element in self.structural_elements:
            try:
                elements.extend(element.faces)
            except AttributeError:
                elements.append(element)
        return elements

    @Attribute
    def primitives(self):
        # Combines mesh subgrids from skin, structure, and loads
        return self.foo + self.structural_element_primitives + self.load_primitives

    @Attribute
    def subgrids(self):
        # Extract unique mesh subgrids from all primitives
        subgrids = (prim.subgrid for prim in self.primitives if hasattr(prim, "subgrid"))
        return list(set(subgrids))  # deduplicates subgrids


# --- Script interface for standalone testing ---
if __name__ == "__main__":
    comp = CodeAster_primitives()
    print(comp.load_primitives)  # Prints the nodes where loads will be applied
