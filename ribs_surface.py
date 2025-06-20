from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from rib_profile import RibProfile
from reference_frame import Frame
import numpy as np


class RibSurface(GeomBase):
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

    rib_spanwise_position = Input(0.7)
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Part
    def rib_frame(self):
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
    def chord_rib(self):
        return self._calculate_chord_at_position(self.rib_spanwise_position)

    @Part
    def rib_curve(self):
        return RibProfile(
            airfoil_name=self.wing_airfoil_root,
            chord=self.chord_rib,
            thickness_factor=self.wing_thickness_factor_root,
            front_spar_position=self.front_spar_position,
            rear_spar_position=self.rear_spar_position,
            position=rotate(
                translate(
                    self.position,
                    "y", self.rib_spanwise_position * self.wing_semi_span,
                    "x",
                    (self.rib_spanwise_position * self.wing_semi_span * np.tan(
                        radians(self.wing_sweep_leading_edge_planform1)))
                    if self.rib_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1
                    else (self.wing_semi_span_planform1 * np.tan(radians(self.wing_sweep_leading_edge_planform1))) + ((
                                                                                                                                  self.rib_spanwise_position * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                        radians(self.wing_sweep_leading_edge_planform2)))

                    # (self.rib_spanwise_position * self.wing_semi_span * tan(
                    #     radians(self.wing_sweep_leading_edge_planform1* (self.wing_semi_span_planform1/(self.rib_spanwise_position * self.wing_semi_span)) + (((self.rib_spanwise_position * self.wing_semi_span) - self.wing_semi_span_planform1)/(self.rib_spanwise_position * self.wing_semi_span))* self.wing_sweep_leading_edge_planform2))
                    # - (self.wing_semi_span_planform1 * tan(radians(self.wing_sweep_leading_edge_planform2 - self.wing_sweep_leading_edge_planform1)))

                ),
                "y",
                radians(
                    self.wing_twist * (self.rib_spanwise_position * self.wing_semi_span) / self.wing_semi_span
                )
            )
        )

    @Part
    def rib_surface(self):
        return Face(island=[self.rib_curve.uppercurve, self.rib_curve.frontspar,  self.rib_curve.lowercurve,  self.rib_curve.rearspar],
                    transparency=0.8
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

    obj = RibSurface(label="Rib")
    display(obj)