from parapy.core import *
import os
from parapy.exchange import STEPWriter
from parapy.gui import display
from parapy.geom import *
from wingbox import Wingbox
from AVL_analysis import WingAVLAnalysis
from visualisation_arrows import LiftArrowArray
import warnings
from FEM_analysis import optimize_plate_thickness
from wing import WingSurface
from meshing import FinalMesh, MeshGenerator
from FEM_analysis import WingFEM, Writer
from find_nodes import CodeAster_primitives
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

DIR = os.path.expanduser("~/Documents")
# DIR = os.path.dirname(__file__)



warnings.filterwarnings("ignore", category=UserWarning)  # Suppress AVL/FEM warnings

# excel_directory = r"C:\Users\nick2\PycharmProjects\KBE\Parameters.xlsm"
excel_directory = r"C:\Users\raane\Documents\Uni\Master\KBE\Year2\Tutorials\Parameters.xlsm"


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
interpolate_airfoil('whitcomb.dat', 'whitcomb_interpolated.dat', factor=25)




class IntegratedWingAnalysis(Base):

    # Wing Parameters
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")
    wing_root_chord = Input(float(pd.read_excel(excel_directory).iloc[0, 1]))
    wing_middle_chord = Input(float(pd.read_excel(excel_directory).iloc[1, 1]))
    wing_tip_chord = Input(float(pd.read_excel(excel_directory).iloc[2, 1]))
    wing_thickness_factor_root = Input(float(pd.read_excel(excel_directory).iloc[4, 1]))
    wing_thickness_factor_middle = Input(float(pd.read_excel(excel_directory).iloc[5, 1]))
    wing_thickness_factor_tip = Input(float(pd.read_excel(excel_directory).iloc[6, 1]))
    wing_semi_span_planform1 = Input(float(pd.read_excel(excel_directory).iloc[8, 1]))
    wing_semi_span = Input(float(pd.read_excel(excel_directory).iloc[9, 1]))
    wing_sweep_leading_edge_planform1 = Input(float(pd.read_excel(excel_directory).iloc[11, 1]))
    wing_sweep_leading_edge_planform2 = Input(float(pd.read_excel(excel_directory).iloc[12, 1]))

    front_spar_position = Input(float(pd.read_excel(excel_directory).iloc[14, 1]))
    rear_spar_position = Input(float(pd.read_excel(excel_directory).iloc[15, 1]))

    rib_number = Input(int(pd.read_excel(excel_directory).iloc[17, 1]))

    stringer_thickness = Input(float(pd.read_excel(excel_directory).iloc[19, 1]))
    stringer_number = Input(int(pd.read_excel(excel_directory).iloc[20, 1]))

    target_deflection = Input(float(pd.read_excel(excel_directory).iloc[23, 5]))
    initial_thickness = Input(float(pd.read_excel(excel_directory).iloc[24, 5]))
    thickness_bounds = Input([float(pd.read_excel(excel_directory).iloc[25, 5]), float(pd.read_excel(excel_directory).iloc[26, 5])])

    # Results Storage
    # avl_results = Attribute()
    # fem_results = Attribute()
    # optimized_parameters = Attribute()
    points_number = Input(int(pd.read_excel(excel_directory).iloc[22, 1]))
    section_number = Input(int(pd.read_excel(excel_directory).iloc[23, 1]))

    Mach = Input(float(pd.read_excel(excel_directory).iloc[17, 5]))
    rho = Input(float(pd.read_excel(excel_directory).iloc[20, 5]))

    element_length = float(pd.read_excel(excel_directory).iloc[27, 5])

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
                           mach=self.Mach,
                           points_number=self.points_number
                           )

    # Stap 2: AVL-analyse (pas zichtbaar in GUI wanneer aangeklikt)
    @Attribute
    def avl_analysis(self):
        return WingAVLAnalysis(
            aircraft=self.wing_surface
            # (label="wing")
            ,
            case_settings=[("alpha_5deg", {'alpha': 5.0})],
            points_number=self.points_number,
            rho=self.rho,
            Mach=self.Mach,
            check_nodes=CodeAster_primitives(wing_airfoil_root=self.wing_airfoil_root,
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
                                             rib_number=self.rib_number,

                                             section_number=self.section_number,
                                             points_number=self.points_number,

                                             finalmesh=self.mesh))

    @Attribute
    def avl_lift_forces_normalized(self):
        lift_forces = []
        max_value = max(self.avl_analysis.lift_forces)
        for i in (self.avl_analysis.lift_forces):
            norm = i / max_value
            lift_forces.append(norm)
        return lift_forces

    @Attribute
    def avl_lift_forces(self):
        return self.avl_analysis.lift_forces

    # @Attribute
    # def torsionbox(self):
    #     return TorsionBox(hidden=True)
    #
    # @Attribute
    # def points_list(self):
    #     return [pt.point for pt in self.torsionbox.points]

    # @Part
    # def lift_arrows(self):
    #     return LiftArrowArray(points_list=self.points_list, lift_forces=self.avl_lift_forces)

    @Attribute
    def points_list(self):
        return self.wing_surface.points

    @Part
    def lift_arrows(self):
        return LiftArrowArray(points_list=[pt.point for pt in self.points_list],
                              lift_forces=self.avl_lift_forces_normalized)

    # FEM Analysis Code

    @Part
    def mesh(self):
        return FinalMesh(check_element=0,
                         wing_airfoil_root=self.wing_airfoil_root,
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
                         rib_number=self.rib_number,

                         section_number=self.section_number,
                         points_number=self.points_number,
                         mesh_generator=MeshGenerator(element_length=self.element_length, shape_to_mesh=self.mesh.shape_to_mesh)
                         )

    @Attribute
    def fem_setup(self):
        return WingFEM(finalmesh=self.mesh,
                       avl=self.avl_analysis,
                       skin_writer=CodeAster_primitives(wing_airfoil_root=self.wing_airfoil_root,
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
                                                        rib_number=self.rib_number,

                                                        section_number=self.section_number,
                                                        points_number=self.points_number,

                                                        finalmesh=self.mesh))

    @Attribute
    def fem_writer(self):
        return Writer(instance=self.fem_setup, avl=self.avl_analysis)

    @Attribute
    def run_fem_analysis(self):
        result = optimize_plate_thickness(
            target_deflection=self.target_deflection,
            wing_fem_instance=self.fem_setup,  # Pass your FEM instance
            writer_instance=self.fem_writer,  # Pass your Writer instance
            initial_thickness=self.initial_thickness,
            thickness_bounds=self.thickness_bounds
        )
        print(f"Optimized thickness: {result['optimized_thickness']} m")
        print(f"Max deflection achieved: {result['max_deflection']} m")
        print(f"Optimization success: {result['success']}")
        return result

    # Showcase optimized thicknesses, and generate stepfile
    @Part
    def wingbox(self):
        return Wingbox(
            wing_airfoil_root=self.wing_airfoil_root,
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
            front_spar_thickness=float(self.run_fem_analysis['optimized_thickness']),
            front_spar_position=self.front_spar_position,
            rear_spar_thickness=float(self.run_fem_analysis['optimized_thickness']),
            rear_spar_position=self.rear_spar_position,

            rib_thickness=float(self.run_fem_analysis['optimized_thickness']),
            rib_number=self.rib_number,

            stringer_thickness=self.stringer_thickness,
            stringer_number=self.stringer_number
        )

    # Step file creation here!!!

    @Attribute
    def check_shape(self):
        shape = self.wingbox
        print("Valid:", shape.shape is not None)
        print("Bounding Box:", shape.bbox if hasattr(shape, 'bbox') else "No bbox")
        return shape

    @Part
    def step_writer_fused(self):
        return STEPWriter(
            default_directory=DIR,
            nodes=[self.check_shape]
        )

if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
    )
    display(obj)
