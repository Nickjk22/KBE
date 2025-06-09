from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np


class SparProfile(Polygon):
    airfoil_name = Input("whitcomb.dat")
    chord = Input(12)
    thickness_factor = Input(1)

    spar_thickness = Input(0.2)
    spar_position = Input(0.2)

    @Attribute
    def points(self):
        with open(self.airfoil_name, 'r') as f:
            points_list = []
            for line in f:
                x, z = line.split(' ', 1)
                if float(self.spar_position) <= float(x) <= float(self.spar_position) + float(self.spar_thickness):
                    points_list.append(self.position.translate(
                        "x", float(x) * self.chord,
                        "z", float(z) * self.chord * self.thickness_factor))
        return points_list


    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

if __name__ == '__main__':
    from parapy.gui import display
    obj = SparProfile(label="Spar profile")
    display(obj)
