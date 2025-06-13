from math import radians, tan
from parapy.geom import *
from parapy.core import *
from reference_frame import Frame
from spars_profile_1d import SparProfile1D
from wing_profiles import WingProfiles


class SparProfile(Wire):
    airfoil_name = Input("whitcomb_interpolated.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    spar_thickness = Input(0.2)
    spar_position = Input(0.2)

    # @Attribute
    # def points(self):
    #     with open(self.airfoil_name, 'r') as f:
    #         points_list = []
    #         for line in f:
    #             x, z = line.split(' ', 1)
    #             if float(self.spar_position) <= float(x) <= float(self.spar_position) + self.spar_thickness/float(self.chord):
    #                 points_list.append(self.position.translate(
    #                     "x", float(x) * self.chord,
    #                     "z", float(z) * self.chord * self.thickness_factor))
    #     return points_list
    #
    #
    # @Part
    # def airfoil_frame(self):
    #     return Frame(pos=self.position,
    #                  hidden=False)
    #



    @Part
    def wing_profiles(self):
        return WingProfiles(airfoil_name=self.airfoil_name,
                            chord=self.chord,
                            thickness_factor=self.thickness_factor,
                            front_spar_position=self.spar_position,
                            rear_spar_position=self.spar_position + self.spar_thickness/self.chord,
                            hidden=True
                            )

    @Part
    def spars_profiles(self):
        return SparProfile1D(airfoil_name=self.airfoil_name,
                             chord=self.chord,
                             thickness_factor=self.thickness_factor,
                             front_spar_position=self.spar_position,
                             rear_spar_position=self.spar_position + self.spar_thickness / self.chord,
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

    @Attribute
    def curves_in(self):
        return [self.uppercurve, self.rearspar, self.lowercurve, self.frontspar]

if __name__ == '__main__':
    from parapy.gui import display
    obj = SparProfile(label="Spar profile")
    display(obj)
