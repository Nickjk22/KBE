from parapy.geom import *
from parapy.core import *
from wing_profiles import WingProfiles
from reference_frame import Frame


class UpperPlateProfile(Wire):
    airfoil_name = Input("whitcomb_interpolated.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    front_spar_thickness = Input(0.1)
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Part
    def wing_profiles(self):
        return WingProfiles(
            airfoil_name=self.airfoil_name,
            chord=self.chord,
            thickness_factor=self.thickness_factor,
            front_spar_position=self.front_spar_position + self.front_spar_thickness / self.chord,
            rear_spar_position=self.rear_spar_position,
            hidden=True
        )

    @Attribute
    def upper_outer_curve(self):
        return self.wing_profiles.upper_curve

    @Attribute
    def upper_inner_curve(self):
        return self.upper_outer_curve.translated("z", -self.front_spar_thickness)

    @Attribute
    def front_upper_connector(self):
        return LineSegment(
            start=self.upper_inner_curve.start,
            end=self.upper_outer_curve.start
        )

    @Attribute
    def rear_upper_connector(self):
        return LineSegment(
            start=self.upper_inner_curve.end,
            end=self.upper_outer_curve.end
        )

    @Attribute
    def upper_plate_wire(self):
        return [self.upper_inner_curve,
                self.rear_upper_connector,
                self.upper_outer_curve,
                self.front_upper_connector]

    @Attribute
    def curves_in(self):
        # Return either upper or lower plate profile here. To visualize both, split into two parts.
        return self.upper_plate_wire  # or self.lower_plate_wire


class LowerPlateProfile(Wire):
    airfoil_name = Input("whitcomb_interpolated.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    front_spar_thickness = Input(0.1)
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Part
    def wing_profiles(self):
        return WingProfiles(
            airfoil_name=self.airfoil_name,
            chord=self.chord,
            thickness_factor=self.thickness_factor,
            front_spar_position=self.front_spar_position + self.front_spar_thickness / self.chord,
            rear_spar_position=self.rear_spar_position,
            hidden=True
        )

    @Attribute
    def lower_outer_curve(self):
        return self.wing_profiles.lower_curve

    @Attribute
    def lower_inner_curve(self):
        return self.lower_outer_curve.translated("z", self.front_spar_thickness)

    @Attribute
    def front_lower_connector(self):
        return LineSegment(
            start=self.lower_inner_curve.start,
            end=self.lower_outer_curve.start
        )

    @Attribute
    def rear_lower_connector(self):
        return LineSegment(
            start=self.lower_inner_curve.end,
            end=self.lower_outer_curve.end
        )

    @Attribute
    def lower_plate_wire(self):
        return [self.lower_inner_curve,
                self.rear_lower_connector,
                self.lower_outer_curve,
                self.front_lower_connector]

    @Attribute
    def curves_in(self):
        # Return either upper or lower plate profile here. To visualize both, split into two parts.
        return self.lower_plate_wire  # or self.lower_plate_wire


if __name__ == '__main__':
    from parapy.gui import display
    obj = LowerPlateProfile(label="Upper Plate Profile")
    display(obj)
