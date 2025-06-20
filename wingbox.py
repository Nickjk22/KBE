from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from wing import WingSurface
from reference_frame import Frame
from spars import Spars
from ribs import Rib
from stringer import Stringer
from sections import Section
import numpy as np
from upper_lower_plates import Plates


class Wingbox(GeomBase):
    # Wing
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

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

    # Spars
    front_spar_thickness = Input(0.1)
    front_spar_position = Input(0.2)
    rear_spar_thickness = Input(0.1)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_thickness = Input(0.1)
    rib_number = Input(15)

    # Stringers
    stringer_thickness = Input(0.01)
    stringer_number = Input(10)

    plate_thickness = Input(0.1)

    # Sections
    section_number = Input(30)

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # Spanwise points for ribs (fraction of semi_span)
    @Attribute
    def spanwise_points_list(self):
        return np.linspace(0, 1, self.rib_number)

    @Attribute
    def spanwise_points_list_sections(self):
        return np.linspace(0, 1, self.section_number)

    # Chordwise points for stringers (fraction of chord)
    @Attribute
    def chordwise_points_list(self):
        return np.linspace(0.1, 0.9, self.stringer_number)

    # Wing airfoil profiles
    @Part
    def wing_root_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_root,
                       chord=self.wing_root_chord,
                       thickness_factor=self.wing_thickness_factor_root,
                       hidden=True)

    @Part
    def wing_middle_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_middle,
                       chord=self.wing_middle_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=translate(self.position, "y", self.wing_semi_span_planform1,
                                          "x", self.wing_semi_span_planform1 * tan(
                               radians(self.wing_sweep_leading_edge_planform1))),
                       hidden=True)

    @Part
    def wing_tip_airfoil(self):
        return Airfoil(airfoil_name=self.wing_airfoil_tip,
                       chord=self.wing_tip_chord,
                       thickness_factor=self.wing_thickness_factor_tip,
                       position=translate(self.position,
                                          "y", self.wing_semi_span,
                                          "x",
                                          self.wing_semi_span_planform1 * np.tan(radians(
                                              self.wing_sweep_leading_edge_planform1)) + (
                                                  (self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                                              radians(
                                                  self.wing_sweep_leading_edge_planform2)))

                                          #                   tan(radians(
                                          # (self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform1 + (1 - self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform2))
                                          ),
                       hidden=True
                       )

    @Part
    def wing_surface(self):
        return WingSurface(wing_airfoil_root=self.wing_airfoil_root,
                           wing_airfoil_middle=self.wing_airfoil_middle,
                           wing_airfoil_tip=self.wing_airfoil_tip,

                           wing_root_chord=self.wing_root_chord,
                           wing_middle_chord=self.wing_middle_chord,
                           wing_tip_chord=self.wing_tip_chord,

                           wing_thickness_factor_root=self.wing_thickness_factor_root,
                           wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                           wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                           wing_semi_span_planform1=self.wing_semi_span_planform1,
                           wing_semi_span=self.wing_semi_span,
                           wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                           wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                           )

    # Spars
    @Part
    def spars(self):
        return Spars(wing_airfoil_root=self.wing_airfoil_root,
                     wing_airfoil_middle=self.wing_airfoil_middle,
                     wing_airfoil_tip=self.wing_airfoil_tip,

                     wing_root_chord=self.wing_root_chord,
                     wing_middle_chord=self.wing_middle_chord,
                     wing_tip_chord=self.wing_tip_chord,

                     wing_thickness_factor_root=self.wing_thickness_factor_root,
                     wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                     wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                     wing_semi_span_planform1=self.wing_semi_span_planform1,
                     wing_semi_span=self.wing_semi_span,
                     wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                     wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                     front_spar_thickness=self.front_spar_thickness,
                     front_spar_position=self.front_spar_position,

                     rear_spar_thickness=self.rear_spar_thickness,
                     rear_spar_position=self.rear_spar_position,
                     )

    # Ribs
    @Part
    def ribs(self):
        return Rib(wing_airfoil_root=self.wing_airfoil_root,
                   wing_airfoil_middle=self.wing_airfoil_middle,
                   wing_airfoil_tip=self.wing_airfoil_tip,

                   wing_root_chord=self.wing_root_chord,
                   wing_middle_chord=self.wing_middle_chord,
                   wing_tip_chord=self.wing_tip_chord,

                   wing_thickness_factor_root=self.wing_thickness_factor_root,
                   wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                   wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                   wing_semi_span_planform1=self.wing_semi_span_planform1,
                   wing_semi_span=self.wing_semi_span,
                   wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                   wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                   rib_thickness=self.rib_thickness,
                   quantify=self.rib_number,
                   rib_spanwise_position=self.spanwise_points_list[child.index],
                   )

    @Part
    def plates(self):
        return Plates(wing_airfoil_root=self.wing_airfoil_root,
                      wing_airfoil_middle=self.wing_airfoil_middle,
                      wing_airfoil_tip=self.wing_airfoil_tip,

                      wing_root_chord=self.wing_root_chord,
                      wing_middle_chord=self.wing_middle_chord,
                      wing_tip_chord=self.wing_tip_chord,

                      wing_thickness_factor_root=self.wing_thickness_factor_root,
                      wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                      wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                      wing_semi_span_planform1=self.wing_semi_span_planform1,
                      wing_semi_span=self.wing_semi_span,
                      wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                      wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                      plate_thickness=self.plate_thickness,
                      front_spar_position=self.front_spar_position,
                      rear_spar_position=self.rear_spar_position,
                      )

    @Part
    def stringers(self):
        return Stringer(wing_airfoil_root=self.wing_airfoil_root,
                        wing_airfoil_middle=self.wing_airfoil_middle,
                        wing_airfoil_tip=self.wing_airfoil_tip,

                        wing_root_chord=self.wing_root_chord,
                        wing_middle_chord=self.wing_middle_chord,
                        wing_tip_chord=self.wing_tip_chord,

                        wing_thickness_factor_root=self.wing_thickness_factor_root,
                        wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                        wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                        wing_semi_span_planform1=self.wing_semi_span_planform1,
                        wing_semi_span=self.wing_semi_span,
                        wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                        wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                        stringer_thickness=self.stringer_thickness,
                        quantify=self.stringer_number,
                        stringer_position=self.chordwise_points_list[child.index],
                        )

    @Part
    def sections(self):
        return Section(wing_airfoil_root=self.wing_airfoil_root,
                       wing_airfoil_middle=self.wing_airfoil_middle,
                       wing_airfoil_tip=self.wing_airfoil_tip,

                       wing_root_chord=self.wing_root_chord,
                       wing_middle_chord=self.wing_middle_chord,
                       wing_tip_chord=self.wing_tip_chord,

                       wing_thickness_factor_root=self.wing_thickness_factor_root,
                       wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                       wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                       wing_semi_span_planform1=self.wing_semi_span_planform1,
                       wing_semi_span=self.wing_semi_span,
                       wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                       wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                       section_number=self.section_number,
                       hidden=True
                       )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Wingbox(label="Wingbox")
    display(obj)
