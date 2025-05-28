from math import radians, tan
from parapy.geom import *
from parapy.core import *
from airfoil import Airfoil
from wing_surfaces import WingSurfaces
from reference_frame import Frame
from spars_surface import SparSurface
from ribs_surface import RibSurface
from kbeutils import avl
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
from sections import Section
import operator


import warnings
from parapy.core.validate import LessThanOrEqualTo, GreaterThan, GreaterThanOrEqualTo, Between, LessThan


# Add these imports at the top of torsionbox.py if not already present
from parapy.core import Attribute, Part
from parapy.geom.occ.brep import BRep


# Interpolate whitcomb airfoil
def interpolate_airfoil(input_file, output_file, factor=10):
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



# class MeshGenerator(Base):
#     shape_to_mesh = Input()
#     element_length = Input(0.2)
#
#     @Attribute
#     def quad_faces(self):
#         faces = [f for f in self.shape_to_mesh.faces if len(f.edges) == 4]
#         lst = []
#
#         taper_ratio = 0.2
#         for f in faces:
#             e1 = f.edges[0]
#             try:
#                 e3 = e1.opposite_edge
#             except Exception:
#                 continue
#
#             if abs(e1.length - e3.length) > taper_ratio * e1.length:
#                 continue
#
#             e2 = f.edges[1]
#             e4 = e2.opposite_edge
#             if abs(e2.length - e4.length) > taper_ratio * e2.length:
#                 continue
#
#             lst.append(f)
#         return lst
#
#     @Part(in_tree=False)
#     def fixed_length(self):
#         return FixedLength(shape_to_mesh=self.shape_to_mesh,
#                           length=self.element_length)
#
#     @Part(in_tree=False)
#     def quad(self):
#         return Quad(quantify=len(self.quad_faces),
#                     shape=self.quad_faces[child.index])
#
#     @Part
#     def tri(self):
#         return Tri(shape_to_mesh=self.shape_to_mesh,
#                   quad_dominant=False,
#                   only_2d=True,
#                   min_size=0.1,
#                   max_size=0.3)
#
#     @Part
#     def mesh(self):
#         return Mesh(shape_to_mesh=self.shape_to_mesh,
#                     display_mode="shaded",
#                     controls=[self.fixed_length, self.quad, self.tri])



Excel_directory = r"C:\Users\raane\Documents\Uni\Master\KBE\Year2\Tutorials\Form.xlsx"
# Input(float(pd.read_excel(Excel_directory).iloc[1, 1]), validator=GreaterThanOrEqualTo(0))

# Class
# Class
class TorsionBox(Base):
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
    front_spar_thickness = Input(0.05)
    front_spar_position = Input(0.2)
    rear_spar_thickness = Input(0.05)
    rear_spar_position = Input(0.6)

    # Ribs
    rib_thickness = Input(0.2)
    rib_number = Input(20)
    section_number = Input(14)

    # Stringers
    # stringer_thickness = Input(0.01)
    # stringer_number = Input(10)

    @Part
    def wing_frame(self):
        return Frame(pos=self.position,
                     hidden=True)

    @Attribute
    def spanwise_points_list(self):
        return np.linspace(0, 1, self.rib_number)

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
                       position=rotate(translate(self.position, "y", self.wing_semi_span_planform1,
                                                 "x", self.wing_semi_span_planform1 * tan(
                               radians(self.wing_sweep_leading_edge_planform1))), "y", radians(
                           self.wing_twist * (self.wing_semi_span_planform1 / self.wing_semi_span))),
                       hidden=True)

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
                                       "y", radians(self.wing_twist)),
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
                            wing_twist=self.wing_twist,

                            front_spar_position=self.front_spar_position,
                            front_spar_thickness=self.front_spar_thickness,

                            rear_spar_position=self.rear_spar_position,
                            rear_spar_thickness=self.rear_spar_thickness,
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
                           wing_twist=self.wing_twist,

                           front_spar_thickness=self.front_spar_thickness,
                           front_spar_position=self.front_spar_position,
                           rear_spar_thickness=self.rear_spar_thickness,
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
                          wing_twist=self.wing_twist,

                          front_spar_position=self.front_spar_position,
                          front_spar_thickness=self.front_spar_thickness,
                          rear_spar_position=self.rear_spar_position,
                          rear_spar_thickness=self.rear_spar_thickness,
                          quantify=self.rib_number,
                          rib_spanwise_position=self.spanwise_points_list[child.index],
                          )

    # # AVL required parts/attributes
    # @Attribute
    # def planform_area(self):  # We assume a trapezoid shape for the wing
    #     area1 = (0.5 * (self.wing_root_chord + self.wing_middle_chord) * self.wing_semi_span_planform1)
    #     area2 = (0.5 * (self.wing_middle_chord + self.wing_tip_chord) * (
    #             self.wing_semi_span - self.wing_semi_span_planform1))
    #     return area1 + area2
    #
    # @Attribute
    # def mac(self):
    #     return (self.wing_semi_span_planform1 / self.wing_semi_span) * (
    #             0.5 * (self.wing_root_chord + self.wing_middle_chord)) + (
    #             (self.wing_semi_span - self.wing_semi_span_planform1) / self.wing_semi_span) * (
    #             0.5 * (self.wing_middle_chord + self.wing_tip_chord))
    #
    # @Attribute
    # def avl_surfaces(self):  # this scans the product tree and collect all instances of the avl.Surface class
    #     return self.find_children(lambda o: isinstance(o, avl.Surface))
    #
    # @Part
    # def avl_configuration(self):
    #     return avl.Configuration(name='aircraft',
    #                              reference_area=self.planform_area,
    #                              reference_span=(self.wing_semi_span * 2),
    #                              reference_chord=self.mac,
    #                              reference_point=self.position.point,
    #                              surfaces=self.avl_surfaces,
    #                              mach=self.mach)


    # @Part
    # def shape_to_mesh(self):
    #     return GeneralFuse(
    #         tools=[
    #                   self.wing_upper_surface,
    #                   self.wing_lower_surface,
    #                   self.front_spar,
    #                   self.rear_spar
    #               ] + [rib.rib_surface for rib in self.ribs],
    #
    #         transparency=0.5
    #     )
    #
    # @Part
    # def mesh_generator(self):
    #     return MeshGenerator(shape_to_mesh=self.shape_to_mesh)
    #
    # @Attribute
    # def mesh(self):
    #     return self.mesh_generator.mesh

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

                       section_number=self.section_number
                       )


if __name__ == '__main__':
    from parapy.gui import display

    interpolate_airfoil('whitcomb.dat', 'whitcomb_interpolated.dat', factor=25)

    obj = TorsionBox(label="Torsion Box")
    display(obj)