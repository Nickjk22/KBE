from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
from spars_profile_1d import SparProfile1D
import numpy as np


class SparSurface(GeomBase):
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

    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Part
    def spars_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # Root airfoil and spars profiles
    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root)

    @Part
    def spar_root_profiles(self):
        return SparProfile1D(airfoil_name=self.wing_airfoil_root,
                             chord=self.wing_root_chord,
                             thickness_factor=self.wing_thickness_factor_root,
                             front_spar_position=self.front_spar_position,
                             rear_spar_position=self.rear_spar_position
                             )

    @Attribute
    def front_spar_root_profile(self):
        return self.spar_root_profiles.front_spar

    @Attribute
    def rear_spar_root_profile(self):
        return self.spar_root_profiles.rear_spar

    # Middle airfoil and spars profiles
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
    def spar_middle_profiles(self):
        return SparProfile1D(airfoil_name=self.wing_airfoil_middle,
                             chord=self.wing_middle_chord,
                             thickness_factor=self.wing_thickness_factor_middle,
                             front_spar_position=self.front_spar_position,
                             rear_spar_position=self.rear_spar_position,
                             position=translate(self.position, "y", self.wing_semi_span_planform1,
                                                "x", self.wing_semi_span_planform1 * tan(
                                     radians(self.wing_sweep_leading_edge_planform1)))
                             )

    @Attribute
    def front_spar_middle_profile(self):
        return self.spar_middle_profiles.front_spar

    @Attribute
    def rear_spar_middle_profile(self):
        return self.spar_middle_profiles.rear_spar

    # Tip airfoil and spars profiles
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
                                          )
                       )

    @Part
    def spar_tip_profiles(self):
        return SparProfile1D(airfoil_name=self.wing_airfoil_tip,
                             chord=self.wing_tip_chord,
                             thickness_factor=self.wing_thickness_factor_tip,
                             front_spar_position=self.front_spar_position,
                             rear_spar_position=self.rear_spar_position,
                             position=translate(self.position,
                                                       "y", self.wing_semi_span,
                                                       "x",
                                                       self.wing_semi_span_planform1 * np.tan(radians(
                                                           self.wing_sweep_leading_edge_planform1)) + (
                                                               (
                                                                       self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                                                           radians(
                                                               self.wing_sweep_leading_edge_planform2)))

                                                       #                   tan(radians(
                                                       # (self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform1 + (1 - self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform2))
                                                       ),
                             )

    @Attribute
    def front_spar_tip_profile(self):
        return self.spar_tip_profiles.front_spar

    @Attribute
    def rear_spar_tip_profile(self):
        return self.spar_tip_profiles.rear_spar

    # Spars surfaces

    @Part
    def front_spar_plan1(self):
        return LoftedSurface(profiles=[self.front_spar_root_profile, self.front_spar_middle_profile],
                             color="Blue",
                             hidden=True)

    @Part
    def rear_spar_plan1(self):
        return LoftedSurface(profiles=[self.rear_spar_root_profile, self.rear_spar_middle_profile],
                             color="Blue",
                             hidden=True)

    @Part
    def front_spar_plan2(self):
        return LoftedSurface(profiles=[self.front_spar_middle_profile, self.front_spar_tip_profile],
                             color="Blue",
                             hidden=True)

    @Part
    def rear_spar_plan2(self):
        return LoftedSurface(profiles=[self.rear_spar_middle_profile, self.rear_spar_tip_profile],
                             color="Blue",
                             hidden=True)


if __name__ == '__main__':
    from parapy.gui import display

    obj = SparSurface(label="Spar")
    display(obj)
