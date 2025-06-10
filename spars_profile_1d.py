from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class SparProfile1D(GeomBase):
    airfoil_name = Input("whitcomb_interpolated.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)


    @Attribute
    def points(self):
        with open(self.airfoil_name, 'r') as f:
            # Read and parse all points
            all_points = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

            # Separate upper and lower surfaces
            upper_points = [(x, z) for x, z in all_points if z >= 0]
            lower_points = [(x, z) for x, z in all_points if z < 0]

            upper_in_region = [(x, z) for x, z in upper_points
                                if self.front_spar_position <= float(x) <= self.rear_spar_position]
            lower_in_region = [(x, z) for x, z in lower_points
                                if self.front_spar_position <= float(x) <= self.rear_spar_position]

            front_spar = []
            rear_spar = []

            # upper coordinate
            front_spar.append(self.position.translate(
                    "x", float(upper_in_region[-1][0]) * self.chord,
                    "z", float(upper_in_region[-1][1]) * self.chord))

            # lower coordinate
            front_spar.append(self.position.translate(
                "x", float(lower_in_region[0][0]) * self.chord,
                "z", float(lower_in_region[0][1]) * self.chord))

            # upper coordinate
            rear_spar.append(self.position.translate(
                "x", float(upper_in_region[0][0]) * self.chord,
                "z", float(upper_in_region[0][1]) * self.chord))

            # lower coordinate
            rear_spar.append(self.position.translate(
                "x", float(lower_in_region[-1][0]) * self.chord,
                "z", float(lower_in_region[-1][1]) * self.chord))

        return {
            'front_spar': front_spar,
            'rear_spar': rear_spar
        }

    @Part
    def front_spar(self):
        return LineSegment(
            start=self.points['front_spar'][0],
            end=self.points['front_spar'][1],
            color="red",
            transparency=0.5
        )

    @Part
    def rear_spar(self):
        return LineSegment(
            start=self.points['rear_spar'][0],
            end=self.points['rear_spar'][1],
            color="red",
            transparency=0.5
        )

    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.airfoil_name,
                       chord=self.chord,
                       thickness_factor=self.thickness_factor)


if __name__ == '__main__':
    from parapy.gui import display

    obj = SparProfile1D(label="Spar profile")
    display(obj)