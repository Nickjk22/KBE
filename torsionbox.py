from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from wing_surfaces import WingSurfaces
from reference_frame import Frame
from spars_surface import SparSurface
from ribs_surface import RibSurface
import numpy as np
from sections import Section
from points import Points
from parapy.core import Attribute, Part



def generate_warning(warning_header, msg):
    from tkinter import Tk, messagebox

    # initialization
    window = Tk()
    window.withdraw()

    # generates message box and waits for user to close it
    messagebox.showwarning(warning_header, msg)

    # kills the gui
    window.deiconify()
    window.destroy()
    window.quit()


# Class
class TorsionBox(Base):
    # Wing
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

    wing_root_chord = Input(6)
    wing_middle_chord = Input(4)
    wing_tip_chord = Input(1.5)

    wing_thickness_factor_root = Input(1)
    wing_thickness_factor_middle = Input(1)
    wing_thickness_factor_tip = Input(1)

    wing_semi_span_planform1 = Input(5)
    wing_semi_span = Input(16)
    wing_sweep_leading_edge_planform1 = Input(20)
    wing_sweep_leading_edge_planform2 = Input(20)

    # Spars
    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_number = Input(12)
    section_number = Input(14)
    points_number = Input(14)

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=True)

    @Attribute
    def spanwise_points_list_ribs(self):
        return np.linspace(0, 1, self.rib_number)

    @Attribute
    def spanwise_points_list(self):
        return np.linspace(0, 1, self.points_number)

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

    # Wing surfaces
    @Part
    def wing_surfaces(self):
        return WingSurfaces(wing_airfoil_root=self.wing_airfoil_root,
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

                            front_spar_position=self.front_spar_position,
                            rear_spar_position=self.rear_spar_position,
                            hidden=True
                            )

    @Part
    def wing_upper_surface(self):
        return SewnShell([self.wing_surfaces.upper_surface1, self.wing_surfaces.upper_surface2])

    @Part
    def wing_lower_surface(self):
        return SewnShell([self.wing_surfaces.lower_surface1, self.wing_surfaces.lower_surface2])

    # Spars
    @Part
    def spars(self):
        return SparSurface(wing_airfoil_root=self.wing_airfoil_root,
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

                           front_spar_position=self.front_spar_position,
                           rear_spar_position=self.rear_spar_position,
                           hidden=True

                           )

    @Part
    def front_spar(self):
        return SewnShell([self.spars.front_spar_plan1, self.spars.front_spar_plan2])

    @Part
    def rear_spar(self):
        return SewnShell([self.spars.rear_spar_plan1, self.spars.rear_spar_plan2])

    # Ribs
    @Part
    def ribs(self):
        return RibSurface(wing_airfoil_root=self.wing_airfoil_root,
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

                          front_spar_position=self.front_spar_position,
                          rear_spar_position=self.rear_spar_position,
                          quantify=self.rib_number,
                          rib_spanwise_position=self.spanwise_points_list_ribs[child.index],
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

    @Part
    def points(self):
        return Points(wing_airfoil_root=self.wing_airfoil_root,
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

                      quantify=self.points_number,
                      point_spanwise_position=self.spanwise_points_list[child.index],
                      )


if __name__ == '__main__':
    from parapy.gui import display
    obj = TorsionBox(label="Torsion Box")
    display(obj)
