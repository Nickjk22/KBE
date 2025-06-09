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
from scipy.interpolate import interp1d


# Interpolate whitcomb airfoil
def interpolate_airfoil(input_file, output_file, factor=5):
    # Read and parse the original data
    with open(input_file, 'r') as f:
        coords = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

    # Split into upper and lower surfaces
    # Find index where y first becomes negative (start of lower surface)
    transition_idx = next(i for i, (x, y) in enumerate(coords) if y < 0)

    upper = coords[:transition_idx]  # Upper surface (including leading edge)
    lower = coords[transition_idx:]  # Lower surface (including trailing edge)

    # Reverse the lower surface for proper parameterization
    lower = lower[::-1]

    # Create separate interpolation functions for upper and lower
    x_upper, y_upper = zip(*upper)
    x_lower, y_lower = zip(*lower)

    # Generate new parameter values (5x denser)
    t_upper = np.linspace(0, 1, len(upper))
    t_lower = np.linspace(0, 1, len(lower))

    new_t = np.linspace(0, 1, len(upper) * factor)

    # Interpolate both x and y coordinates
    fx_upper = interp1d(t_upper, x_upper, kind='cubic')
    fy_upper = interp1d(t_upper, y_upper, kind='cubic')

    fx_lower = interp1d(t_lower, x_lower, kind='cubic')
    fy_lower = interp1d(t_lower, y_lower, kind='cubic')

    # Generate new coordinates
    new_upper = list(zip(fx_upper(new_t), fy_upper(new_t)))
    new_lower = list(zip(fx_lower(new_t), fy_lower(new_t)))[::-1]  # Reverse back

    # Combine while maintaining correct order (upper TE -> LE -> lower LE -> TE)
    new_coords = new_upper + new_lower[1:]

    # Write to file (maintaining original format)
    with open(output_file, 'w') as f:
        for x, y in new_coords:
            f.write(f"{x:.5f} {y:.5f}\n")

    print(f"Created {output_file} with {len(new_coords)} points")


# Usage
interpolate_airfoil('whitcomb.dat', 'whitcomb_interpolated.dat', factor=10)


# Class
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
    wing_twist = Input(0)

    # Spars
    front_spar_thickness = Input(1)
    front_spar_position = Input(0.2)
    rear_spar_thickness = Input(1)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_thickness = Input(0.2)
    rib_number = Input(15)

    # Stringers
    stringer_thickness = Input(0.01)
    stringer_number = Input(10)

    # Sections
    section_number = Input(30)

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=False)

    # @Attribute
    # def spanwise_points_list_plan1(self):
    #     return np.linspace(0, 1, self.rib_number_plan1)
    #
    # @Attribute
    # def spanwise_points_list_plan2(self):
    #     return np.linspace(0, 1, self.rib_number_plan2)

    # Spanwise points for ribs (fraction of semi_span)
    @Attribute
    def spanwise_points_list(self):
        return np.linspace(0, 1, self.rib_number)

    @Attribute
    def spanwise_points_list_sections(self):
        return np.linspace(0, 1, self.section_number)

    #Chordwise points for stringers (fraction of chord)
    @Attribute
    def chordwise_points_list(self):
        return np.linspace(0.1, 0.9, self.stringer_number)

    # Wing airfoil profiles
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
                                                         (self.wing_semi_span - self.wing_semi_span_planform1) * np.tan(
                                                     radians(
                                                         self.wing_sweep_leading_edge_planform2)))

                                                 #                   tan(radians(
                                                 # (self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform1 + (1 - self.wing_semi_span_planform1/self.wing_semi_span)*self.wing_sweep_leading_edge_planform2))
                                                 ),
                                       "y", radians(self.wing_twist))
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
                           wing_twist=self.wing_twist

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
                     wing_twist=self.wing_twist,

                     front_spar_thickness=self.front_spar_thickness,
                     front_spar_position=self.front_spar_position,

                     rear_spar_thickness=self.rear_spar_thickness,
                     rear_spar_position=self.rear_spar_position,
                     )

    # @Part
    # def ribs_plan1(self):
    #     return Rib(
    #         wing_airfoil_root=self.wing_airfoil_root,
    #         wing_airfoil_tip=self.wing_airfoil_middle,
    #
    #         wing_root_chord=self.wing_root_chord,
    #         wing_tip_chord=self.wing_middle_chord,
    #
    #         wing_thickness_factor_root=self.wing_thickness_factor_root,
    #         wing_thickness_factor_tip=self.wing_thickness_factor_middle,
    #
    #         wing_semi_span=self.wing_semi_span_planform1,
    #         wing_sweep=self.wing_sweep_leading_edge_planform1,
    #         wing_twist=self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span),
    #
    #         rib_thickness=self.rib_thickness,
    #         quantify=self.rib_number_plan1,
    #         rib_spanwise_position=self.spanwise_points_list_plan1[child.index]
    #     )
    #
    # @Part
    # def ribs_plan2(self):
    #     return Rib(position=rotate(translate(self.position, "y", self.wing_semi_span_planform1,
    #                                          "x", self.wing_semi_span_planform1 * tan(
    #             radians(self.wing_sweep_leading_edge_planform1))), "y", radians(
    #         self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span))),
    #                wing_airfoil_root=self.wing_airfoil_middle,
    #                wing_airfoil_tip=self.wing_airfoil_tip,
    #
    #                wing_root_chord=self.wing_middle_chord,
    #                wing_tip_chord=self.wing_tip_chord,
    #
    #                wing_thickness_factor_root=self.wing_thickness_factor_middle,
    #                wing_thickness_factor_tip=self.wing_thickness_factor_tip,
    #
    #                wing_semi_span=self.wing_semi_span - self.wing_semi_span_planform1,
    #                wing_sweep=(self.wing_sweep_leading_edge_planform2 + self.wing_sweep_leading_edge_planform1) / 2,
    #                wing_twist=self.wing_twist * (1 - self.wing_semi_span_planform1 / self.wing_semi_span),
    #
    #                rib_thickness=self.rib_thickness,
    #                quantify=self.rib_number_plan2,
    #                rib_spanwise_position=self.spanwise_points_list_plan2[child.index]
    #                )

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
                   wing_twist=self.wing_twist,

                   rib_thickness=self.rib_thickness,
                   quantify=self.rib_number,
                   rib_spanwise_position=self.spanwise_points_list[child.index],
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
                        wing_twist=self.wing_twist,

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
                       wing_twist=self.wing_twist,

                       section_number=self.section_number,
                       hidden=True
                       )


if __name__ == '__main__':
    from parapy.gui import display

    obj = Wingbox(label="Wingbox")
    display(obj)
