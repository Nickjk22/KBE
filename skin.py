from typing import List, Sequence

from parapy.core import Attribute, Base, Input, Part
from parapy.core import *
from aircraft_fem.examples.aircraft.thickness import DEFAULT_SKIN_THICKNESS, THICKNESSES, \
    THICKNESS_TO_RGB
from parapy.core.widgets import Dropdown
from aircraft_fem.examples.aircraft.material import STEEL
from parapy.gui.events import EVT_SELECTION_CHANGING
from parapy.geom import FittedCurve, IntersectedShapes, LoftedSolid, \
    ModifiedShape, Plane, Point, RuledSurface, VZ
from parapy.core.decorators import on_event
from torsionbox import TorsionBox
from parapy.geom import LineSegment
from parapy.geom import Point
from sections import Section
from meshing import FinalMesh


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
    section_number = Input(13)
    wing_semi_span = Input(30)

    @Attribute
    def torsionbox(self):
        return TorsionBox()

    @Part
    def sections(self):
        return Section()

    @Attribute
    def cp_points(self):
        positions = []
        for tellen in range(self.section_number):
            x_tel = self.sections.section_airfoil[tellen].chord * 1 / 4 + \
                    self.sections.section_airfoil[tellen].position[0]
            y_tel = (self.wing_semi_span / self.section_number) * tellen

            point = Point(x_tel, y_tel, 0)

            positions.append(point)

        return positions

    @Part
    def finalmesh(self):
        return FinalMesh()

    @Part
    def skin(self):
        return Skin(faces_to_keep=self.torsionbox.wing_upper_surface,
                    wing_shell=self.finalmesh.shape_to_mesh)

    @Part
    def segments(self):
        return LineSegment(
            start=self.cp_points[child.index].translate("z", -100),
            end=self.cp_points[child.index].translate("z", 100),
            quantify=len(self.cp_points),
            color="green"  # optional for visual clarity
        )



    @Attribute
    def strip_force_meshable_shape_intersection_pts(self):
        pts = []

        points = self.cp_points

        for i, point in enumerate(points):
            segment = LineSegment(point.translate("z", -100), point.translate("z", 100))

            for face in self.skin.shell.faces:
                try:
                    pt = face.intersection_point(segment)
                    pts.append(pt)
                    break
                except RuntimeError:
                    continue
            else:
                # Fallback at start and end
                if i == 0:
                    pts.append(point)
                elif i == len(points) - 1:
                    pts.append(point)
                else:
                    raise RuntimeError(f"No intersection at {segment}, revise input or geometry.")
        return pts


if __name__ == '__main__':
    from parapy.gui import display

    obj = CodeAster_primitives()

    display(obj)
