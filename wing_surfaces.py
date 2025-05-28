from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np
from wing_profiles import WingProfiles


class WingSurfaces(GeomBase):
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

    front_spar_position = Input(0.2)
    front_spar_thickness = Input(0.1)

    rear_spar_position = Input(0.6)
    rear_spar_thickness = Input(0.1)

    @Attribute
    def profiles_upper1(self):
        return [self.wing_upper_root, self.wing_upper_middle]

    @Attribute
    def profiles_upper2(self):
        return [self.wing_upper_middle, self.wing_upper_tip]

    @Attribute
    def profiles_lower1(self):
        return [self.wing_lower_root, self.wing_lower_middle]

    @Attribute
    def profiles_lower2(self):
        return [self.wing_lower_middle, self.wing_lower_tip]

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # Root
    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root)

    @Part
    def wing_root(self):
        return WingProfiles(airfoil_name=self.wing_airfoil_root,
                            chord=self.wing_root_chord,
                            thickness_factor=self.wing_thickness_factor_root,

                            front_spar_position=self.front_spar_position,
                            front_spar_thickness=self.front_spar_thickness,

                            rear_spar_position=self.rear_spar_position,
                            rear_spar_thickness=self.rear_spar_thickness
                            )

    @Attribute
    def wing_upper_root(self):
        return self.wing_root.upper_curve

    @Attribute
    def wing_lower_root(self):
        return self.wing_root.lower_curve

    # Middle
    @Part
    def wing_middle_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_middle,
                       chord=self.wing_middle_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       hidden=True,
                       position=rotate(translate(self.position, "y", self.wing_semi_span_planform1,
                                                 "x", self.wing_semi_span_planform1 * tan(
                               radians(self.wing_sweep_leading_edge_planform1))), "y", radians(
                           self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span))))

    @Part
    def wing_middle(self):
        return WingProfiles(airfoil_name=self.wing_airfoil_middle,
                            chord=self.wing_middle_chord,
                            thickness_factor=self.wing_thickness_factor_middle,

                            front_spar_position=self.front_spar_position,
                            front_spar_thickness=self.front_spar_thickness,

                            rear_spar_position=self.rear_spar_position,
                            rear_spar_thickness=self.rear_spar_thickness,
                            position=rotate(translate(self.position, "y", self.wing_semi_span_planform1,
                                                      "x", self.wing_semi_span_planform1 * tan(
                                    radians(self.wing_sweep_leading_edge_planform1))), "y", radians(
                                self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span)))
                            )

    @Attribute
    def wing_upper_middle(self):
        return self.wing_middle.upper_curve

    @Attribute
    def wing_lower_middle(self):
        return self.wing_middle.lower_curve

    # Tip
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

    @Part
    def wing_tip(self):
        return WingProfiles(airfoil_name=self.wing_airfoil_tip,
                            chord=self.wing_tip_chord,
                            thickness_factor=self.wing_thickness_factor_tip,

                            front_spar_position=self.front_spar_position,
                            front_spar_thickness=self.front_spar_thickness,

                            rear_spar_position=self.rear_spar_position,
                            rear_spar_thickness=self.rear_spar_thickness,
                            position=rotate(translate(self.position,
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
                                            "y", radians(self.wing_twist))
                            )

    @Attribute
    def wing_upper_tip(self):
        return self.wing_tip.upper_curve

    @Attribute
    def wing_lower_tip(self):
        return self.wing_tip.lower_curve

    # Surfaces
    @Part
    def upper_surface1(self):
        return LoftedSurface(profiles=self.profiles_upper1,
                             color="Blue",
                             transparency=1,
                             hidden=True)

    @Part
    def upper_surface2(self):
        return LoftedSurface(profiles=self.profiles_upper2,
                             color="Blue",
                             transparency=1,
                             hidden=True)

    @Part
    def lower_surface1(self):
        return LoftedSurface(profiles=self.profiles_lower1,
                             color="Blue",
                             transparency=1,
                             hidden=True)

    @Part
    def lower_surface2(self):
        return LoftedSurface(profiles=self.profiles_lower2,
                             color="Blue",
                             transparency=1,
                             hidden=True)
    #
    # @Part
    # def wing_upper_surface(self):
    #     return FusedShell(
    #         shape_in=self.upper_surface1,
    #         tool=[self.upper_surface2],
    #         mesh_deflection=0.0001,
    #         transparency=0.8,
    #         color="Yellow"
    #     )
    #
    # @Part
    # def wing_lower_surface(self):
    #     return FusedShell(
    #         shape_in=self.lower_surface1,
    #         tool=[self.lower_surface2],
    #         mesh_deflection=0.0001,
    #         transparency=0.8,
    #         color="Yellow"
    #     )

if __name__ == '__main__':
    from parapy.gui import display

    obj = WingSurfaces(label="wing surface")
    display(obj)
