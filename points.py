from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class Points(GeomBase):
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

    point_spanwise_position = Input(0.7)

    @Part
    def point_frame(self):
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

    @Attribute
    def y_at_x_025(self):
        # Load airfoil data
        data = np.loadtxt(self.wing_airfoil_root)

        number = int((len(data[:, 0]) + 1) / 2)
        # slice rows 0 … mid-1  (upper surface)
        x_coords = data[:number, 0]  # first-half x column
        y_coords = data[:number, 1]  # first-half y column

        # Find index of x closest to 0.25
        index_closest = np.argmin(np.abs(x_coords - 0.25))

        # Now we want to find all y-values with *exactly* the same x
        x_target = x_coords[index_closest]
        matching_y = y_coords[np.isclose(x_coords, x_target)]

        # Return the largest of them
        return max(matching_y)

    @Attribute
    def chord_point(self):
        return self._calculate_chord_at_position(self.point_spanwise_position)

    @Attribute
    def point(self):
        return Point(((self.point_spanwise_position * self.wing_semi_span * np.tan(
            radians(self.wing_sweep_leading_edge_planform1)))
                      if self.point_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1
                      else (self.wing_semi_span_planform1 * np.tan(
            radians(self.wing_sweep_leading_edge_planform1))) + ((
                                                                         self.point_spanwise_position * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
            radians(self.wing_sweep_leading_edge_planform2)))) + self.chord_point * 0.25,

                     self.point_spanwise_position * self.wing_semi_span,

                     (self.chord_point * self.y_at_x_025))

    @Part
    def visualized_point(self):
        return Vertex(
            point=Point(
                ((self.point_spanwise_position * self.wing_semi_span * np.tan(
                    radians(self.wing_sweep_leading_edge_planform1)))
                 if self.point_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1
                 else (self.wing_semi_span_planform1 * np.tan(
                    radians(self.wing_sweep_leading_edge_planform1))) + (
                              (
                                      self.point_spanwise_position * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                          radians(self.wing_sweep_leading_edge_planform2))))
                + self.chord_point * 0.25,

                self.point_spanwise_position * self.wing_semi_span,

                (self.chord_point * self.y_at_x_025)
            ),
            label="Quarter Chord Point"
        )

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root,
                       hidden=True)

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
                                       "y", radians(self.wing_twist)),
                       hidden=True
                       )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Points(label="Points")
    display(obj)
