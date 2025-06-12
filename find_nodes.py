from typing import List, Sequence

from parapy.core import Attribute, Base, Input, Part
from parapy.core import *
from thickness import DEFAULT_SKIN_THICKNESS, THICKNESSES, \
    THICKNESS_TO_RGB
from parapy.core.widgets import Dropdown
from material import STEEL
from parapy.gui.events import EVT_SELECTION_CHANGING
from parapy.geom import FittedCurve, IntersectedShapes, LoftedSolid, \
    ModifiedShape, Plane, Point, RuledSurface, VZ
from parapy.core.decorators import on_event
from torsionbox import TorsionBox
from parapy.geom import LineSegment
from parapy.geom import Point
from sections import Section
from meshing import MeshGenerator
from points import Points
import numpy as np





class Skin(Base):
    wing_shell = Input()
    faces_to_keep = Input()
    thickness = Input(DEFAULT_SKIN_THICKNESS, widget=Dropdown(THICKNESSES))
    material = Input(STEEL)

    @on_event(EVT_SELECTION_CHANGING)
    def handle_selection_clicks(self, evt):
        """I am the owner of all clicks on my subshapes"""
        evt.owner = self

    @Part
    def shell(self):
        return ModifiedShape(self.wing_shell,
                             keep=self.faces_to_keep,
                             color=THICKNESS_TO_RGB[self.thickness])


class CodeAster_primitives(Base):
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

    # Spars
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_number = Input(12)

    # Extra
    section_number = Input(14)
    points_number = Input(14)

    element_length = Input(0.1)

    @Input
    def torsionbox(self):
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
                          rear_spar_position =self.rear_spar_position,
                          rib_number = self.rib_number,

                          section_number=self.section_number,
                          points_number = self.points_number,
                          hidden=True)

    @Input
    def finalmesh(self):
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
        return Skin(faces_to_keep=self.torsionbox.wing_upper_surface,
                    wing_shell=self.finalmesh.shape_to_mesh)

    @Attribute
    def spanwise_points_list(self):
        return np.linspace(0, 1, self.points_number)

    # @Part
    # def points(self):
    #     return Points(wing_airfoil_root=self.wing_airfoil_root,
    #                   wing_airfoil_middle=self.wing_airfoil_middle,
    #                   wing_airfoil_tip=self.wing_airfoil_tip,
    #
    #                   wing_root_chord=self.wing_root_chord,
    #                   wing_middle_chord=self.wing_middle_chord,
    #                   wing_tip_chord=self.wing_tip_chord,
    #
    #                   wing_thickness_factor_root=self.wing_thickness_factor_root,
    #                   wing_thickness_factor_middle=self.wing_thickness_factor_middle,
    #                   wing_thickness_factor_tip=self.wing_thickness_factor_tip,
    #
    #                   wing_semi_span_planform1=self.wing_semi_span_planform1,
    #                   wing_semi_span=self.wing_semi_span,
    #                   wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
    #                   wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,
    #                   wing_twist=self.wing_twist,
    #
    #                   quantify=self.points_number,
    #                   point_spanwise_position=self.spanwise_points_list[child.index],
    #                   hidden=False
    #                   )

    @Attribute
    def points_list(self):
        return [pt.point for pt in self.torsionbox.points]

    @Attribute
    def load_primitives(self):
        lst = []
        mesh = self.finalmesh.mesh
        for pt in self.points_list:
            label = 'group_no_' + str(id(pt))
            # node = mesh.grid.find_node_at(pt)
            tolerance = self.finalmesh.element_length * 2
            node = mesh.grid.find_node_at(pt,
                                          tolerance=tolerance)
            if node is None:
                msg = "No mesh node found near strip force cp point. Increase" \
                      "the tolerance of find_node_at."
                raise RuntimeError(msg)

            # subgrid = NodesSubGrid(nodes=[node],
            #                        label=label)
            # load = NodeLoad(subgrid.nodes,
            #                 subgrid=subgrid,
            #                 outward=True,
            #                 force=Vector(0, 0, round(strip_result.lift_force, 2)),
            #                 base_length=strip_result.arrow.base_length,
            #                 head_length=strip_result.arrow.head_length)
            lst.append(node)
        return lst

    @Attribute
    def foo(self):
        lst = []
        history = self.finalmesh.shape_to_mesh.history
        mesh = self.finalmesh.mesh
        get_subgrid = mesh.get_subgrid_on_the_fly

        for face in self.skin.shell.faces:
            shape = history(face)
            label = 'group_ma_' + str(id(face))
            subgrid = get_subgrid(shape[0], *shape[1:], label=label)
            lst.append(subgrid)

        return lst

    # @Attribute
    # def foo(self):
    #     lst = []
    #     history = self.finalmesh.shape_to_mesh.history
    #     mesh = self.finalmesh.mesh_generator.mesh
    #     get_subgrid = mesh.get_subgrid_on_the_fly
    #
    #     for face in self.root_skin.shell.faces:
    #         shape = history(face)
    #         label = 'group_ma_' + str(id(face))
    #         subgrid = get_subgrid(shape[0], *shape[1:], label=label)
    #         lst.append(subgrid)
    #
    #     for skin in [self.tip_skin, self.lateral_skin]:
    #         for face in skin.shell.faces:
    #             shapes = history(face)
    #             label = 'group_ma_' + str(id(face))
    #             if len(shapes) == 1:
    #                 subgrid = get_subgrid(shapes[0], label=label)
    #             else:
    #                 subgrid = get_subgrid(shapes[0], *shapes[1:], label=label)
    #             lst.append(subgrid)
    #     return lst

    @Attribute
    def structural_elements(self):
        return (
                self.torsionbox.ribs +
                [self.torsionbox.front_spar] +
                [self.torsionbox.rear_spar] +
                [self.torsionbox.wing_upper_surface] +
                [self.torsionbox.wing_lower_surface]
        )

    @property
    def structural_element_primitives(self):
        elements = []
        for element in self.structural_elements:
            try:
                elements.extend(element.faces)
            except AttributeError:
                elements.append(element)
        return elements

    @Attribute
    def primitives(self):
        return self.foo + self.structural_element_primitives + self.load_primitives

    @Attribute
    def subgrids(self):
        subgrids = (prim.subgrid for prim in self.primitives if hasattr(prim, "subgrid"))
        return list(set(subgrids))  # removes duplicates, e.g. root clamps and root shell use same subgrids


if __name__ == "__main__":
    comp = CodeAster_primitives()

    print(comp.load_primitives)

