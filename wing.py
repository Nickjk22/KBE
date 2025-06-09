from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from reference_frame import Frame
import numpy as np
import kbeutils.avl as avl



class WingSurface(GeomBase):
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
    wing_twist = Input(-10)

    mach = Input(0.4)

    @Attribute
    def profiles(self):
        return [self.wing_root_airfoil, self.wing_middle_airfoil]

    @Attribute
    def profiles2(self):
        return [self.wing_middle_airfoil, self.wing_tip_airfoil]

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root)

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
    def wing_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.wing_tip_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=rotate(translate(self.position,
                                                 "y", self.wing_semi_span,
                                                 "x",
                                                 self.wing_semi_span_planform1 * np.tan(radians(
                                                     self.wing_sweep_leading_edge_planform1)) + (
                                                         (self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(radians(
                                                     self.wing_sweep_leading_edge_planform2)))

                                                 #                   tan(radians(
                                                 # (self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform1 + (1 - self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform2))
                                                 ),
                                       "y", radians(self.wing_twist))
                       )

    @Part
    def lofted_solid(self):
        return LoftedSurface(profiles=self.profiles,
                             color="Blue",
                             transparency=1,
                             hidden=True)

    @Part
    def lofted_solid2(self):
        return LoftedSurface(profiles=self.profiles2,
                             color="Blue",
                             transparency=1,
                             hidden=True)

    @Part
    def right_wing(self):
        return FusedShell(
            shape_in=self.lofted_solid2,
            tool=[self.lofted_solid],
            mesh_deflection=0.0001,
            transparency=0.8,
            color="White"
        )

    @Attribute
    def planform_area(self):
        """Berekening van het planform-oppervlak van een samengestelde vleugel (trapeziumvormig)"""
        span1 = self.wing_semi_span_planform1
        span2 = self.wing_semi_span - span1

        area1 = 0.5 * (self.wing_root_chord + self.wing_middle_chord) * span1
        area2 = 0.5 * (self.wing_middle_chord + self.wing_tip_chord) * span2

        return area1 + area2

    @Attribute
    def mac(self):
        """Berekening van de Mean Aerodynamic Chord (MAC) over twee vleugeldelen"""
        span1 = self.wing_semi_span_planform1
        span2 = self.wing_semi_span - span1

        mac1 = 0.5 * (self.wing_root_chord + self.wing_middle_chord)
        mac2 = 0.5 * (self.wing_middle_chord + self.wing_tip_chord)

        return (span1 * mac1 + span2 * mac2) / self.wing_semi_span

    # AVL
    @Attribute
    def wing_sections(self):
        return [self.wing_root_airfoil, self.wing_middle_airfoil, self.wing_tip_airfoil]

    @Part
    def avl_surface(self):
        return avl.Surface(name=self.name,
                           n_chordwise=12,
                           chord_spacing=avl.Spacing.cosine,
                           n_spanwise=20,
                           span_spacing=avl.Spacing.cosine,
                           y_duplicate=None,
                           transparency=1,
                           sections=[Airfoil.avl_section
                                     for Airfoil in self.wing_sections])

    @Attribute
    def avl_surfaces(self):  # this scans the product tree and collect all instances of the avl.Surface class
        return self.find_children(lambda o: isinstance(o, avl.Surface))

    @Part
    def avl_configuration(self):
        return avl.Configuration(name='wing',
                                 reference_area=self.planform_area,
                                 reference_span=self.wing_semi_span,
                                 reference_chord=self.mac,
                                 reference_point=self.position.point,
                                 surfaces=self.avl_surfaces,
                                 mach=self.mach)


if __name__ == '__main__':
    from parapy.gui import display

    obj = WingSurface(label="wing surface")
    display(obj)
