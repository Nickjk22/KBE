from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np
from spars_profile_1d import SparProfile1D
from wing_profiles import WingProfiles


class RibProfile(GeomBase):
    airfoil_name = Input("whitcomb.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    front_spar_position = Input(0.2)
    front_spar_thickness = Input(0.1)

    rear_spar_position = Input(0.6)
    rear_spar_thickness = Input(0.1)

    @Attribute
    def points(self):
        with open(self.airfoil_name, 'r') as f:
            points_list = []
            for line in f:
                x, z = line.split(' ', 1)
                if self.front_spar_position + 0.5 * self.front_spar_thickness <= float(
                        x) <= self.rear_spar_position + 0.5 * self.rear_spar_thickness:
                    points_list.append(self.position.translate(
                        "x", float(x) * self.chord,
                        "z", float(z) * self.chord * self.thickness_factor))
        return points_list

    @Part
    def wing_profiles(self):
        return WingProfiles(airfoil_name=self.airfoil_name,
                            chord=self.chord,
                            thickness_factor=self.thickness_factor,

                            front_spar_position=self.front_spar_position,
                            front_spar_thickness=self.front_spar_thickness,

                            rear_spar_position=self.rear_spar_position,
                            rear_spar_thickness=self.rear_spar_thickness,
                            hidden=True
                            )

    @Part
    def spars_profiles(self):
        return SparProfile1D(airfoil_name=self.airfoil_name,
                             chord=self.chord,
                             thickness_factor=self.thickness_factor,
                             front_spar_thickness=self.front_spar_thickness,
                             front_spar_position=self.front_spar_position,
                             rear_spar_thickness=self.rear_spar_thickness,
                             rear_spar_position=self.rear_spar_position,
                             hidden=True
                             )

    @Attribute
    def uppercurve(self):
        return self.wing_profiles.upper_curve

    @Attribute
    def lowercurve(self):
        return self.wing_profiles.lower_curve

    @Attribute
    def frontspar(self):
        return self.spars_profiles.front_spar

    @Attribute
    def rearspar(self):
        return self.spars_profiles.rear_spar


    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=True)


if __name__ == '__main__':
    from parapy.gui import display

    obj = RibProfile(label="rib profile")
    display(obj)