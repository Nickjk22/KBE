from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
from spar_profile import SparProfile
import numpy as np
from stringer_profile import StringerProfile
from upper_lower_plate_profile import UpperLowerPlateProfile


class Plates(GeomBase):
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

    plate_thickness = Input(0.1)
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Part
    def plate_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # Root airfoil and spars profiles
    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root)

    @Part
    def plates_root_profile(self):
        return UpperLowerPlateProfile(airfoil_name=self.wing_airfoil_root,
                                      chord=self.wing_root_chord,
                                      thickness_factor=self.wing_thickness_factor_root,
                                      front_spar_thickness=self.plate_thickness,
                                      front_spar_position=self.front_spar_position,
                                      rear_spar_position=self.rear_spar_position,
                                      hidden=True)

    # Middle airfoil and stringer profiles
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
    def plates_middle_profile(self):
        return UpperLowerPlateProfile(airfoil_name=self.wing_airfoil_middle,
                                      chord=self.wing_middle_chord,
                                      thickness_factor=self.wing_thickness_factor_middle,
                                      front_spar_thickness=self.plate_thickness,
                                      front_spar_position=self.front_spar_position,
                                      rear_spar_position=self.rear_spar_position,
                                      position=translate(self.position, "y", self.wing_semi_span_planform1,
                                                         "x", self.wing_semi_span_planform1 * tan(
                                              radians(self.wing_sweep_leading_edge_planform1))),
                                      hidden=True
                                      )

    # Tip airfoil and stringer profiles
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
    def plates_tip_profile(self):
        return UpperLowerPlateProfile(airfoil_name=self.wing_airfoil_tip,
                                      chord=self.wing_tip_chord,
                                      thickness_factor=self.wing_thickness_factor_tip,
                                      front_spar_thickness=self.plate_thickness,
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
                                      hidden=True
                                      )

    @Part
    def upper_plate_loft(self):
        return LoftedSolid(
            profiles=[
                self.plates_root_profile.upper_plate,
                self.plates_middle_profile.upper_plate,
            ],
            color="red",
            hidden=True,
        )

    @Part
    def upper_plate_loft2(self):
        return LoftedSolid(
            profiles=[
                self.plates_middle_profile.upper_plate,
                self.plates_tip_profile.upper_plate,
            ],
            color="red",
            hidden=True
        )

    @Part
    def lower_plate_loft(self):
        return LoftedSolid(
            profiles=[
                self.plates_root_profile.lower_plate,
                self.plates_middle_profile.lower_plate,
            ],
            color="red",
            hidden=True
        )

    @Part
    def lower_plate_loft2(self):
        return LoftedSolid(
            profiles=[
                self.plates_middle_profile.lower_plate,
                self.plates_tip_profile.lower_plate,
            ],
            color="red",
            hidden=True
        )

    @Part
    def upper_plate(self):
        return FusedSolid(
            shape_in=self.upper_plate_loft,
            tool=[self.upper_plate_loft2],
            mesh_deflection=0.0001,
            color="Blue",
        )

    @Part
    def lower_plate(self):
        return FusedSolid(
            shape_in=self.lower_plate_loft,
            tool=[self.lower_plate_loft2],
            mesh_deflection=0.0001,
            color="Blue",
        )

    # @Part
    # def lower_stringer_loft(self):
    #     return LoftedSolid(
    #         profiles=[
    #             self.stringer_root_profile.lower_stringer,
    #             self.stringer_middle_profile.lower_stringer,
    #             self.stringer_tip_profile.lower_stringer
    #         ],
    #         color="blue",
    #         transparency=0.5
    #     )

    # @Attribute
    # def profiles_plan1(self):
    #     return [
    #         self.stringer_root_profile.upper_stringer,
    #         self.stringer_middle_profile.upper_stringer
    #     ]
    #
    # @Attribute
    # def profiles_plan2(self):
    #     return [
    #         self.stringer_middle_profile.upper_stringer,
    #         self.stringer_tip_profile.upper_stringer
    #     ]

    # @Part
    # def stringer_plan1(self):
    #     return LoftedSolid(
    #         profiles=self.profiles_plan1,
    #         color="Blue",
    #         transparency=0.5
    #     )
    #
    # @Part
    # def stringer_plan2(self):
    #     return LoftedSolid(
    #         profiles=self.profiles_plan2,
    #         color="Blue",
    #         transparency=0.5
    #     )
    #
    # @Part
    # def stringer(self):
    #     return FusedSolid(
    #         shape_in=self.stringer_plan1,
    #         tool=[self.stringer_plan2],
    #         color="Blue",
    #         transparency=0.5
    #     )

    # @Part
    # def upper_stringer_loft(self):
    #     return LoftedSolid(
    #         profiles=[
    #             self.stringer_root_profile.upper_stringer,
    #             self.stringer_middle_profile.upper_stringer,
    #             self.stringer_tip_profile.upper_stringer
    #         ],
    #         color="red",
    #         transparency=0.3
    #     )
    #
    # @Part
    # def lower_stringer_loft(self):
    #     return LoftedSolid(
    #         profiles=[
    #             self.stringer_root_profile.lower_stringer,
    #             self.stringer_middle_profile.lower_stringer,
    #             self.stringer_tip_profile.lower_stringer
    #         ],
    #         color="blue",
    #         transparency=0.3
    #     )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Plates(label="Plates")
    display(obj)
