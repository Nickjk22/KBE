from parapy.geom import *
from parapy.core import *
from reference_frame import Frame
import kbeutils.avl as avl


class Airfoil(FittedCurve):
    airfoil_name = Input()
    chord = Input(1)
    thickness_factor = Input(1)
    mesh_deflection = 0.0001

    @Attribute
    def points(self):
        with open(self.airfoil_name, 'r') as f:
            points_list = []
            for line in f:
                x, z = line.split(' ', 1)
                points_list.append(self.position.translate(
                    "x", float(x) * self.chord,
                    "z", float(z) * self.chord * self.thickness_factor))
        return points_list

    @Part
    def airfoil_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # AVL required parts/attributes
    @Part
    def avl_section(self):
        return avl.SectionFromCurve(curve_in=self, hidden=True)


if __name__ == '__main__':
    from parapy.gui import display
    obj = Airfoil(label="airfoil", airfoil_name="whitcomb.dat")
    display(obj)
