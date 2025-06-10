from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame


class Rib(LoftedSolid):
    wing_airfoil_root = Input("whitcomb.dat")
    wing_airfoil_tip = Input("whitcomb.dat")

    wing_root_chord = Input(12)
    wing_tip_chord = Input(3)

    wing_thickness_factor_root = Input(1)
    wing_thickness_factor_tip = Input(1)

    wing_sweep = Input(20)
    wing_semi_span = Input(30)

    rib_thickness = Input(0.2)
    rib_spanwise_position = Input(0.5)

    @Attribute
    def profiles(self):
        return [self.rib_root_airfoil, self.rib_tip_airfoil]

    @Part
    def rib_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    @Attribute
    def root_chord_rib(self):
        return self.rib_spanwise_position * self.wing_semi_span * (
                (self.wing_tip_chord - self.wing_root_chord) / self.wing_semi_span) + self.wing_root_chord

    @Attribute
    def tip_chord_rib(self):
        return (self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness) * (
                (self.wing_tip_chord - self.wing_root_chord) / self.wing_semi_span) + self.wing_root_chord

    @Part
    def rib_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.root_chord_rib,
                       thickness_factor=self.wing_thickness_factor_root,
                       position=translate(self.position,
                                          "y", self.rib_spanwise_position * self.wing_semi_span,
                                          "x", self.rib_spanwise_position * self.wing_semi_span * tan(radians(
                               self.wing_sweep)))
                       )

    @Part
    def rib_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.tip_chord_rib,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=translate(self.position,
                                          "y",
                                          self.rib_spanwise_position * self.wing_semi_span + self.rib_thickness,
                                          "x", (
                                                  self.wing_semi_span * self.rib_spanwise_position + self.rib_thickness) * tan(
                               radians(
                                   self.wing_sweep)))
                       )

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root,
                       position=self.position,
                       hidden=True
                       )

    @Part
    def wing_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.wing_tip_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=translate(self.position, "y", self.wing_semi_span,
                                          "x", self.wing_semi_span * tan(
                               radians(self.wing_sweep)),
                                          ),
                       hidden=True)

    @Part
    def lofted_solid(self):
        return LoftedSolid(profiles=self.profiles,
                           color="Blue",
                           hidden=not (__name__ == '__main__'))


if __name__ == '__main__':
    from parapy.gui import display

    obj = Rib(label="Rib", mesh_deflection=0.0001)
    display(obj)
