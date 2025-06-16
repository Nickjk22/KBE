from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame


class WingProfiles(GeomBase):
    airfoil_name = Input("whitcomb.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    @Attribute
    def wing_surface_profiles(self):
        with open(self.airfoil_name, 'r') as f:
            # Read and parse all points
            all_points = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

            # Separate upper and lower surfaces
            upper_points = [(x, z) for x, z in all_points if z >= 0]
            lower_points = [(x, z) for x, z in all_points if z < 0]

            # Get points within stringer region
            upper_in_region = [(x, z) for x, z in upper_points
                               if self.front_spar_position <= x <= self.rear_spar_position]
            lower_in_region = [(x, z) for x, z in lower_points
                               if self.front_spar_position <= x <= self.rear_spar_position]

            # Sort points by x-coordinate
            upper_in_region.sort(key=lambda p: p[0])
            lower_in_region.sort(key=lambda p: p[0])

            # Create upper stringer profile (rectangle)
            upper_profile = []
            for x, z in upper_in_region:
                upper_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", z * self.chord * self.thickness_factor))

            # Create lower stringer profile (rectangle)
            lower_profile = []
            for x, z in lower_in_region:
                lower_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", z * self.chord * self.thickness_factor))

        return {
            'upper_profile': upper_profile,
            'lower_profile': lower_profile
        }

    @Part
    def airfoil(self):
        return Airfoil(airfoil_name=self.airfoil_name,
                       chord=self.chord,
                       thickness_factor=self.thickness_factor,
                       hidden=True)

    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    @Part
    def upper_curve(self):
        return FittedCurve(
            points=self.wing_surface_profiles['upper_profile'],
            color="red",
            transparency=0.5
        )

    @Part
    def lower_curve(self):
        return FittedCurve(
            points=self.wing_surface_profiles['lower_profile'],
            color="blue",
            transparency=0.5
        )


if __name__ == '__main__':
    from parapy.gui import display
    obj = WingProfiles(label="Wing profiles")
    display(obj)
