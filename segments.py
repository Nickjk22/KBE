from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class Segment(GeomBase):
    wing_airfoil_root = Input("whitcomb.dat")
    wing_airfoil_middle = Input("whitcomb.dat")
    wing_airfoil_tip = Input("whitcomb.dat")

    wing_root_chord = Input(12)
    wing_middle_chord = Input(7)
    wing_tip_chord = Input(3)

    wing_thickness_factor_root = Input(1)
    wing_thickness_factor_middle = Input(1)
    wing_thickness_factor_tip = Input(1)

    wing_semi_span_planform1 = Input(10)
    wing_semi_span = Input(30)
    wing_sweep_leading_edge_planform1 = Input(20)
    wing_sweep_leading_edge_planform2 = Input(20)
    wing_twist = Input(0)

    segment_number = Input(14)

    @Attribute
    def spanwise_points_list_segments(self):
        return np.linspace(0, 1, self.segment_number)

    @Part
    def segment_frame(self):
        return Frame(pos=self.position,
                     hidden=True)

    def _calculate_chord_at_position(self, normalized_position):
        if normalized_position * self.wing_semi_span < self.wing_semi_span_planform1:
            return (normalized_position / (
                    self.wing_semi_span_planform1 / self.wing_semi_span)) * self.wing_semi_span_planform1 * (
                    (
                            self.wing_middle_chord - self.wing_root_chord) / self.wing_semi_span_planform1) + self.wing_root_chord
        else:
            return ((normalized_position - (self.wing_semi_span_planform1 / self.wing_semi_span)) / (
                    1 - (self.wing_semi_span_planform1 / self.wing_semi_span))) * (
                    self.wing_semi_span - self.wing_semi_span_planform1) * (
                    (self.wing_tip_chord - self.wing_middle_chord) / (
                    self.wing_semi_span - self.wing_semi_span_planform1)) + self.wing_middle_chord

    # @Attribute
    # def chord_section(self):
    #     return self._calculate_chord_at_position(self.segment_spanwise_position)

    @Part
    def section_airfoil(self):
        return LineSegment(quantify=self.segment_number,
            start=Point(((self.spanwise_points_list_segments[child.index]
                                   * self.wing_semi_span * np.tan(
                    radians(self.wing_sweep_leading_edge_planform1)))
                                  if self.spanwise_points_list_segments[child.index]
                                     * self.wing_semi_span < self.wing_semi_span_planform1
                                  else (self.wing_semi_span_planform1 * np.tan(
            radians(self.wing_sweep_leading_edge_planform1))) + ((
                                                                         self.spanwise_points_list_segments[child.index]
                                                                         * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
            radians(self.wing_sweep_leading_edge_planform2)))) + self._calculate_chord_at_position(
            self.spanwise_points_list_segments[child.index])
                                 * 0.25, self.spanwise_points_list_segments[child.index]
                                 * self.wing_semi_span,
                                 10),
                           end=Point(((self.spanwise_points_list_segments[child.index]
                                              * self.wing_semi_span * np.tan(
                    radians(self.wing_sweep_leading_edge_planform1)))
                                             if self.spanwise_points_list_segments[child.index]
                                                * self.wing_semi_span < self.wing_semi_span_planform1
                                             else (self.wing_semi_span_planform1 * np.tan(
            radians(self.wing_sweep_leading_edge_planform1))) + ((
                                                                         self.spanwise_points_list_segments[
                                                                             child.index]
                                                                         * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
            radians(self.wing_sweep_leading_edge_planform2)))) + self._calculate_chord_at_position(
            self.spanwise_points_list_segments[child.index])
                                            * 0.25, self.spanwise_points_list_segments[child.index]
                                            * self.wing_semi_span,
                                            -10),

                           )

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root)

    @Part
    def wing_middle_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_middle,
                       chord=self.wing_middle_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=rotate(translate(self.position, "y", self.wing_semi_span_planform1,
                                                 "x", self.wing_semi_span_planform1 * tan(
                               radians(self.wing_sweep_leading_edge_planform1))), "y", radians(
                           self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span))),
                       hidden=True)

    @Part
    def wing_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.wing_tip_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=rotate(translate(self.position,
                                                 "y", self.wing_semi_span,
                                                 "x",
                                                 self.wing_semi_span_planform1 * np.tan(radians(
                                                     self.wing_sweep_leading_edge_planform1)) + (
                                                         (self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                                                     radians(
                                                         self.wing_sweep_leading_edge_planform2)))

                                                 #                   tan(radians(
                                                 # (self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform1 + (1 - self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform2))
                                                 ),
                                       "y", radians(self.wing_twist))
                       )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Segment(label="Segment")
    display(obj)
