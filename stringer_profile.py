from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class StringerProfile(Polygon):
    airfoil_name = Input("whitcomb.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    stringer_thickness = Input(0.05)
    stringer_position = Input(0.2)

    @Attribute
    def stringer_profiles(self):
        with open(self.airfoil_name, 'r') as f:
            # Read and parse all points
            all_points = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

            # Separate upper and lower surfaces
            upper_points = [(x, z) for x, z in all_points if z >= 0]
            lower_points = [(x, z) for x, z in all_points if z < 0]

            # Get points within stringer region
            upper_in_region = [(x, z) for x, z in upper_points
                               if self.stringer_position <= x <= self.stringer_position + self.stringer_thickness]
            lower_in_region = [(x, z) for x, z in lower_points
                               if self.stringer_position <= x <= self.stringer_position + self.stringer_thickness]

            # Sort points by x-coordinate
            upper_in_region.sort(key=lambda p: p[0])
            lower_in_region.sort(key=lambda p: p[0])

            # Create upper stringer profile (rectangle)
            upper_profile = []
            # Forward pass: original upper points
            for x, z in upper_in_region:
                upper_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", z * self.chord * self.thickness_factor))
            # Backward pass: offset points (z - thickness)
            for x, z in reversed(upper_in_region):
                upper_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", (z - self.stringer_thickness/self.thickness_factor) * self.chord * self.thickness_factor))

            # Create lower stringer profile (rectangle)
            lower_profile = []
            # Forward pass: original lower points
            for x, z in lower_in_region:
                lower_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", z * self.chord * self.thickness_factor))
            # Backward pass: offset points (z + thickness)
            for x, z in reversed(lower_in_region):
                lower_profile.append(self.position.translate(
                    "x", x * self.chord,
                    "z", (z + self.stringer_thickness/self.thickness_factor) * self.chord * self.thickness_factor))

        return {
            'upper_profile': upper_profile,
            'lower_profile': lower_profile
        }

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.airfoil_name,
                       chord=self.chord,
                       thickness_factor=self.thickness_factor,
                       hidden=True)
    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    @Part
    def upper_stringer(self):
        return Polygon(
            points=self.stringer_profiles['upper_profile'],
            color="red",
            transparency=0.5
        )

    @Part
    def lower_stringer(self):
        return Polygon(
            points=self.stringer_profiles['lower_profile'],
            color="blue",
            transparency=0.5
        )

if __name__ == '__main__':
    from parapy.gui import display
    obj = StringerProfile(label="Stringer profile")
    display(obj)
