from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class Rib(LoftedSolid):
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

    rib_thickness = Input(0.2)
    rib_spanwise_position = Input(0.7)

    @Attribute
    def profiles(self):
        return [self.rib_root_airfoil, self.rib_tip_airfoil]

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
    def root_chord_rib(self):
        return self._calculate_chord_at_position(self.rib_spanwise_position)

    @Attribute
    def tip_chord_rib(self):
        offset_position = self.rib_spanwise_position + (self.rib_thickness / self.wing_semi_span)
        if offset_position > 1.0:
            offset_position = 1.0
        return self._calculate_chord_at_position(offset_position)

    @Part
    def rib_root_airfoil(self):
        return Airfoil(
            airfoil_name=self.wing_airfoil_root,
            chord=self.root_chord_rib,
            thickness_factor=self.wing_thickness_factor_root,
            position=translate(
                    self.position,
                    "y", self.rib_spanwise_position * self.wing_semi_span,
                    "x",
                    (self.rib_spanwise_position * self.wing_semi_span * np.tan(
                        radians(self.wing_sweep_leading_edge_planform1)))
                    if self.rib_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1
                    else (self.wing_semi_span_planform1 * np.tan(radians(self.wing_sweep_leading_edge_planform1))) + ((
                                                                                                                                  self.rib_spanwise_position * self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                        radians(self.wing_sweep_leading_edge_planform2)))

                    # (self.rib_spanwise_position * self.wing_semi_span * tan( radians(
                # self.wing_sweep_leading_edge_planform1* (self.wing_semi_span_planform1/(self.rib_spanwise_position
                # * self.wing_semi_span)) + (((self.rib_spanwise_position * self.wing_semi_span) -
                # self.wing_semi_span_planform1)/(self.rib_spanwise_position * self.wing_semi_span))*
                # self.wing_sweep_leading_edge_planform2)) - (self.wing_semi_span_planform1 * tan(radians(
                # self.wing_sweep_leading_edge_planform2 - self.wing_sweep_leading_edge_planform1)))
                ),
            hidden=True

            )

    @Part
    def rib_tip_airfoil(self):
        return Airfoil(
            airfoil_name=self.wing_airfoil_root,
            chord=self.tip_chord_rib,
            # if self.rib_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1 else self.root_chord_rib,
            thickness_factor=self.wing_thickness_factor_root,
            position=translate(
                    self.position,
                    "y", self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness,
                    "x",
                    (self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness) * tan(
                        radians(self.wing_sweep_leading_edge_planform1)
                    ) if self.rib_spanwise_position * self.wing_semi_span < self.wing_semi_span_planform1
                    else (self.wing_semi_span_planform1 * tan(radians(self.wing_sweep_leading_edge_planform1))) + (((
                                                                                                                                self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness) - self.wing_semi_span_planform1) * tan(
                        radians(self.wing_sweep_leading_edge_planform2)))

                    # (
                    #         (self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness) * tan(
                    #     radians(self.wing_sweep_leading_edge_planform2))
                    # - self.wing_semi_span_planform1 * tan(radians(self.wing_sweep_leading_edge_planform1))

                ),
            hidden=True
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
                       position=translate(self.position, "y", self.wing_semi_span_planform1,
                                                 "x", self.wing_semi_span_planform1 * tan(
                               radians(self.wing_sweep_leading_edge_planform1))),
                       hidden=True)

    @Part
    def wing_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.wing_tip_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=translate(self.position,
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
                       hidden=True
                       )

    @Part
    def lofted_solid(self):
        return LoftedSolid(profiles=self.profiles,
                           color=[211,211,211],
                          )

    @property
    def volume(self):
        return (self.lofted_solid.volume
        )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Rib(label="Rib", mesh_deflection=0.0001)
    display(obj)
