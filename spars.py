from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
from spar_profile import SparProfile
import numpy as np


class Spars(GeomBase):
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

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

    front_spar_thickness = Input(1)
    front_spar_position = Input(0.2)

    rear_spar_thickness = Input(1)
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
                       thickness_factor=self.wing_thickness_factor_root,
                       hidden=True)

    @Part
    def front_spar_root_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_root,
                           chord=self.wing_root_chord,
                           thickness_factor=self.wing_thickness_factor_root,
                           spar_thickness=self.front_spar_thickness,
                           spar_position=self.front_spar_position,
                           hidden=True)

    @Part
    def rear_spar_root_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_root,
                           chord=self.wing_root_chord,
                           thickness_factor=self.wing_thickness_factor_root,
                           spar_thickness=self.rear_spar_thickness,
                           spar_position=self.rear_spar_position,
                           hidden=True)

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
    def front_spar_middle_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_middle,
                           chord=self.wing_middle_chord,
                           thickness_factor=self.wing_thickness_factor_middle,
                           spar_thickness=self.front_spar_thickness,
                           spar_position=self.front_spar_position,
                           position=translate(self.position, "y", self.wing_semi_span_planform1,
                                              "x", self.wing_semi_span_planform1 * tan(
                                   radians(self.wing_sweep_leading_edge_planform1))),
                           hidden=True
                           )

    @Part
    def rear_spar_middle_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_middle,
                           chord=self.wing_middle_chord,
                           thickness_factor=self.wing_thickness_factor_middle,
                           spar_thickness=self.rear_spar_thickness,
                           spar_position=self.rear_spar_position,
                           position=translate(self.position, "y", self.wing_semi_span_planform1,
                                              "x", self.wing_semi_span_planform1 * tan(
                                   radians(self.wing_sweep_leading_edge_planform1))),
                           hidden=True
                           )

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
                                          ),
                       hidden=True
                       )

    @Part
    def front_spar_tip_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_tip,
                           chord=self.wing_tip_chord,
                           thickness_factor=self.wing_thickness_factor_tip,
                           spar_thickness=self.front_spar_thickness,
                           spar_position=self.front_spar_position,
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
    def rear_spar_tip_profile(self):
        return SparProfile(airfoil_name=self.wing_airfoil_tip,
                           chord=self.wing_tip_chord,
                           thickness_factor=self.wing_thickness_factor_tip,
                           spar_thickness=self.rear_spar_thickness,
                           spar_position=self.rear_spar_position,
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
                           hidden=True
                           )

    @Attribute
    def profiles_front_plan1(self):
        return [self.front_spar_root_profile, self.front_spar_middle_profile]

    @Attribute
    def profiles_front_plan2(self):
        return [self.front_spar_middle_profile, self.front_spar_tip_profile]

    @Attribute
    def profiles_rear_plan1(self):
        return [self.rear_spar_root_profile, self.rear_spar_middle_profile]

    @Attribute
    def profiles_rear_plan2(self):
        return [self.rear_spar_middle_profile, self.rear_spar_tip_profile]

    @Part
    def front_spar_plan1(self):
        return LoftedSolid(profiles=self.profiles_front_plan1,
                           color="Blue",
                           hidden=True)

    @Part
    def front_spar_plan2(self):
        return LoftedSolid(profiles=self.profiles_front_plan2,
                           color="Blue",
                           hidden=True)

    @Part
    def rear_spar_plan1(self):
        return LoftedSolid(profiles=self.profiles_rear_plan1,
                           color="Blue",
                           hidden=True)

    @Part
    def rear_spar_plan2(self):
        return LoftedSolid(profiles=self.profiles_rear_plan2,
                           color="Blue",
                           hidden=True)

    @Part
    def front_spar(self):
        return FusedSolid(
            shape_in=self.front_spar_plan1,
            tool=[self.front_spar_plan2],
            mesh_deflection=0.0001,
            color=[169,169,169],
        )

    @Part
    def rear_spar(self):
        return FusedSolid(
            shape_in=self.rear_spar_plan1,
            tool=[self.rear_spar_plan2],
            mesh_deflection=0.0001,
            color=[169,169,169],
        )

    # @Part
    # def solid_spar(self):
    #     return FusedSolid(
    #         shape_in=Solid(self.front_spar),
    #         tool=[Solid(self.rear_spar)],
    #         color="Cyan",
    #         hidden=True,
    #     )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Spars(label="Spars")
    display(obj)
