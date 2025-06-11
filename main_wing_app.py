from parapy.core import *
import os
from parapy.exchange import STEPWriter
from parapy.gui import display
from parapy.geom import *
from wingbox import Wingbox
from parapy.geom import Compound
from AVL_analysis import WingAVLAnalysis
from visualisation_arrows import LiftArrowArray
import warnings
from FEM_analysis import optimize_plate_thickness
from wing import WingSurface
from meshing import MeshGenerator
from FEM_analysis import WingFEM, Writer
from find_nodes import CodeAster_primitives
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from kbeutils import avl
from parapy.lib.code_aster import (_F, DEFI_MATERIAU)
from parapy.core.validate import LessThanOrEqualTo, GreaterThan, GreaterThanOrEqualTo, Between, LessThan



DIR = os.path.expanduser("~/Documents")
# DIR = os.path.dirname(__file__)


warnings.filterwarnings("ignore", category=UserWarning)  # Suppress AVL/FEM warnings

# excel_directory = r"C:\Users\nick2\PycharmProjects\KBE\Parameters.xlsm"


excel_directory = r"C:\Users\raane\Documents\Uni\Master\KBE\Year2\Tutorials\Parameters.xlsm"

def interpolate_airfoil(input_file, output_file, factor=5):
    # Read and parse the original data
    with open(input_file, 'r') as f:
        coords = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

    # Find index where y first becomes negative (start of lower surface)
    transition_idx = next(i for i, (x, y) in enumerate(coords) if y < 0)

    # Upper and lower surfaces
    upper = coords[:transition_idx]
    lower = coords[transition_idx:]
    lower = lower[::-1]  # Reverse for correct direction (TE to LE)

    x_upper, y_upper = zip(*upper)
    x_lower, y_lower = zip(*lower)

    # Create a common set of x-coordinates from LE to TE (monotonic)
    x_min = max(min(x_upper), min(x_lower))
    x_max = min(max(x_upper), max(x_lower))

    # Interpolate based on a common x-grid
    num_points = max(len(x_upper), len(x_lower)) * factor
    common_x = np.linspace(x_min, x_max, num_points)

    # Interpolate y as a function of x for both surfaces
    f_upper = interp1d(x_upper, y_upper, kind='cubic', fill_value="extrapolate")
    f_lower = interp1d(x_lower, y_lower, kind='cubic', fill_value="extrapolate")

    new_upper = list(zip(common_x, f_upper(common_x)))
    new_lower = list(zip(common_x, f_lower(common_x)))

    # Combine: upper from TE to LE, lower from LE to TE
    new_coords = new_upper[::-1] + new_lower[1:]

    # Write to output file
    with open(output_file, 'w') as f:
        for x, y in new_coords:
            f.write(f"{x:.5f} {y:.5f}\n")

    print(f"Created {output_file} with {len(new_coords)} points and shared x-coordinates")


# Usage
interpolate_airfoil('whitcomb.dat', 'whitcomb_interpolated.dat', factor=10)

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

class IntegratedWingAnalysis(Base):
    # Wing Parameters
    wing_airfoil_root = Input("whitcomb_interpolated.dat")
    wing_airfoil_middle = Input("whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("whitcomb_interpolated.dat")

    # Warnings
    popup_gui_semi_span = Input(True)
    popup_gui_front_spar_position = Input(True)
    popup_gui_rear_spar_position = Input(True)
    popup_gui_rib_number = Input(True)
    popup_gui_thickness_bounds = Input(True)
    popup_gui_initial_thickness = Input(True)
    popup_gui_root_chord = Input(True)
    popup_gui_middle_chord = Input(True)

    wing_root_chord = Input(float(pd.read_excel(excel_directory).iloc[0, 1]), validator=GreaterThan(0))
    wing_middle_chord = Input(float(pd.read_excel(excel_directory).iloc[1, 1]), validator=GreaterThan(0))
    wing_tip_chord = Input(float(pd.read_excel(excel_directory).iloc[2, 1]), validator=GreaterThan(0))

    @Attribute
    def corrected_root_chord(self):
        if self.wing_root_chord < self.wing_middle_chord:
            msg = f"Root chord ({self.wing_root_chord}) cannot be smaller than middle chord ({self.wing_middle_chord}). Input will be ignored and root chord will be set equal to middle chord."
            warnings.warn(msg)
            if self.popup_gui_root_chord:
                generate_warning("Warning: Value changed", msg)
            return self.wing_middle_chord
        else:
            return self.wing_root_chord

    @Attribute
    def corrected_middle_chord(self):
        if self.wing_middle_chord < self.wing_tip_chord:
            msg = f"Middle chord ({self.wing_middle_chord}) cannot be smaller than tip chord ({self.wing_tip_chord}). Input will be ignored and middle chord will be set equal to tip chord."
            warnings.warn(msg)
            if self.popup_gui_middle_chord:
                generate_warning("Warning: Value changed", msg)
            return self.wing_tip_chord
        else:
            return self.wing_middle_chord

    wing_thickness_factor_root = Input(float(pd.read_excel(excel_directory).iloc[4, 1]), validator=GreaterThan(0))
    wing_thickness_factor_middle = Input(float(pd.read_excel(excel_directory).iloc[5, 1]), validator=GreaterThan(0))
    wing_thickness_factor_tip = Input(float(pd.read_excel(excel_directory).iloc[6, 1]), validator=GreaterThan(0))

    wing_semi_span_planform1 = Input(float(pd.read_excel(excel_directory).iloc[8, 1]), validator=GreaterThan(0))
    wing_semi_span = Input(float(pd.read_excel(excel_directory).iloc[9, 1]), validator=GreaterThan(0))

    @Attribute
    def corrected_semi_span(self):
        if self.wing_semi_span < self.wing_semi_span_planform1:
            msg = f"Wing semi-span ({self.wing_semi_span}) should be larger than or equal to the semi-span of planform 1 ({self.wing_semi_span}). Input will be ignored and will be set equal to each other."
            warnings.warn(msg)
            if self.popup_gui_semi_span:
                generate_warning("Warning: Value changed", msg)
            return self.wing_semi_span_planform1
        else:
            return self.wing_semi_span

    wing_sweep_leading_edge_planform1 = Input(float(pd.read_excel(excel_directory).iloc[11, 1]), validator=Between(-60, 60))
    wing_sweep_leading_edge_planform2 = Input(float(pd.read_excel(excel_directory).iloc[12, 1]), validator=Between(-60,60))

    front_spar_position = Input(float(pd.read_excel(excel_directory).iloc[14, 1]), validator=Between(0,1))
    rear_spar_position = Input(float(pd.read_excel(excel_directory).iloc[15, 1]), validator=Between(0,1))

    @Attribute
    def corrected_front_spar_position(self):
        if self.front_spar_position + self.corrected_thickness_bounds[1]/self.wing_tip_chord > self.rear_spar_position:
            msg = f"Front spar position ({self.front_spar_position}) plus its maximum thickness should be smaller than the rear spar position ({self.rear_spar_position}). Input will be ignored and front spar position will be placed to the front."
            warnings.warn(msg)
            if self.popup_gui_front_spar_position:
                generate_warning("Warning: Value changed", msg)
            return self.rear_spar_position-self.corrected_thickness_bounds[1]/self.wing_tip_chord
        else:
            return self.front_spar_position

    @Attribute
    def corrected_rear_spar_position(self):
        if self.rear_spar_position + self.corrected_thickness_bounds[1]/self.wing_tip_chord > 1:
            msg = f"Rear spar position ({self.rear_spar_position}) plus its maximum thickness should be smaller than the chord length. Input will be ignored and rear spar position will be placed to the front."
            warnings.warn(msg)
            if self.popup_gui_rear_spar_position:
                generate_warning("Warning: Value changed", msg)
            return 1-self.corrected_thickness_bounds[1]/self.wing_tip_chord
        else:
            return self.rear_spar_position

    stringer_thickness = Input(float(pd.read_excel(excel_directory).iloc[19, 1]), validator=GreaterThan(0))
    stringer_number = Input(int(pd.read_excel(excel_directory).iloc[20, 1]), validator=GreaterThan(0))

    target_deflection = Input(float(pd.read_excel(excel_directory).iloc[23, 5]), validator=GreaterThan(0))
    thickness_bounds = Input(
        [float(pd.read_excel(excel_directory).iloc[25, 5]), float(pd.read_excel(excel_directory).iloc[26, 5])])
    initial_thickness = Input(float(pd.read_excel(excel_directory).iloc[24, 5]), validator=GreaterThan(0))

    @Attribute
    def corrected_thickness_bounds(self):
        if self.thickness_bounds[0] > self.thickness_bounds[1]:
            msg = f"Lower thickness bound ({self.thickness_bounds[0]}) should be smaller than the upper thickness bound ({self.thickness_bounds[1]}). Input will be ignored and lower bound will be set equal to upper bound."
            warnings.warn(msg)
            if self.popup_gui_thickness_bounds:
                generate_warning("Warning: Value changed", msg)
            return [self.thickness_bounds[1], self.thickness_bounds[1]]
        else:
            return [self.thickness_bounds[0], self.thickness_bounds[1]]

    @Attribute
    def corrected_initial_thickness(self):
        if not self.corrected_thickness_bounds[0] < self.initial_thickness < self.corrected_thickness_bounds[1]:
            msg = f"Initial thickness ({self.initial_thickness}) should be in between lower and upper bound ({self.corrected_thickness_bounds}). Input will be ignored and set to average of bounds."
            warnings.warn(msg)
            if self.popup_gui_initial_thickness:
                generate_warning("Warning: Value changed", msg)
            return (self.corrected_thickness_bounds[0]+self.corrected_thickness_bounds[1])/2
        else:
            return self.initial_thickness


    rib_number = Input(int(pd.read_excel(excel_directory).iloc[17, 1]), validator=GreaterThan(0))

    @Attribute
    def corrected_rib_number(self):
        if self.rib_number*self.thickness_bounds[1] > self.corrected_semi_span:
            msg = f"Number of ribs ({self.rib_number}) times maximum rib thickness (thickness upper bound) ({self.thickness_bounds[1]}) exceeds the wing semi-span. Input will be ignored and the rib number will be set equal to wing semi-span divided by the maximum thickness."
            warnings.warn(msg)
            if self.popup_gui_rib_number:
                generate_warning("Warning: Value changed", msg)
            return self.wing_semi_span/self.thickness_bounds[1]
        else:
            return self.rib_number



    # Results Storage
    # avl_results = Attribute()
    # fem_results = Attribute()
    # optimized_parameters = Attribute()
    points_number = Input(int(pd.read_excel(excel_directory).iloc[22, 1]), validator=GreaterThan(0))
    section_number = Input(int(pd.read_excel(excel_directory).iloc[23, 1]), validator=GreaterThan(0))

    Mach = Input(float(pd.read_excel(excel_directory).iloc[17, 5]), validator=Between(0,1))
    rho = Input(float(pd.read_excel(excel_directory).iloc[20, 5]), validator=GreaterThan(0))

    element_length = Input(float(pd.read_excel(excel_directory).iloc[27, 5]), validator=Between(0,1))
    is_mirrored = Input(True)

    material = Input('Aluminium')

    @Part
    def wing_surface(self):
        return WingSurface(wing_airfoil_root=self.wing_airfoil_root,
                           wing_airfoil_middle=self.wing_airfoil_middle,
                           wing_airfoil_tip=self.wing_airfoil_tip,

                           wing_root_chord=self.corrected_root_chord,
                           wing_middle_chord=self.corrected_middle_chord,
                           wing_tip_chord=self.wing_tip_chord,

                           wing_thickness_factor_root=self.wing_thickness_factor_root,
                           wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                           wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                           wing_semi_span_planform1=self.wing_semi_span_planform1,
                           wing_semi_span=self.corrected_semi_span,
                           wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                           wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,
                           mach=self.Mach,
                           points_number=self.points_number,
                           is_mirrored=self.is_mirrored
                           )

    # Read the settings based case number (1 = fixed CL, 2 = fixed angle of attack)
    new_setting = Input(2)
    new_parameter = Input(5)

    # Write the cases based on the input
    @Attribute
    def case_inputs(self):
        if self.new_setting == 2:
            case1 = ('fixed_aoa', {'alpha': self.new_parameter})
            return case1
        else:
            case2 = ('fixed_cl', {'alpha': avl.Parameter(name='alpha', value=self.new_parameter, setting='CL')})
            return case2

    # Stap 2: AVL-analyse (pas zichtbaar in GUI wanneer aangeklikt)
    @Attribute
    def avl_analysis(self):
        return WingAVLAnalysis(
            aircraft=self.wing_surface
            # (label="wing")
            ,
            case_settings=[self.case_inputs],
            is_mirrored=self.is_mirrored,
            points_number=(2*self.points_number) if self.is_mirrored else self.points_number,
            rho=self.rho,
            Mach=self.Mach,
            check_nodes=CodeAster_primitives(wing_airfoil_root=self.wing_airfoil_root,
                                             wing_airfoil_middle=self.wing_airfoil_middle,
                                             wing_airfoil_tip=self.wing_airfoil_tip,

                                             wing_root_chord=self.corrected_root_chord,
                                             wing_middle_chord=self.corrected_middle_chord,
                                             wing_tip_chord=self.wing_tip_chord,

                                             wing_thickness_factor_root=self.wing_thickness_factor_root,
                                             wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                                             wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                                             wing_semi_span_planform1=self.wing_semi_span_planform1,
                                             wing_semi_span=self.corrected_semi_span,
                                             wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                                             wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                                             front_spar_position=self.corrected_front_spar_position,
                                             rear_spar_position=self.corrected_rear_spar_position,
                                             rib_number=self.corrected_rib_number,

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
        return MeshGenerator(check_element=0,
                             wing_airfoil_root=self.wing_airfoil_root,
                             wing_airfoil_middle=self.wing_airfoil_middle,
                             wing_airfoil_tip=self.wing_airfoil_tip,

                             wing_root_chord=self.corrected_root_chord,
                             wing_middle_chord=self.corrected_middle_chord,
                             wing_tip_chord=self.wing_tip_chord,

                             wing_thickness_factor_root=self.wing_thickness_factor_root,
                             wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                             wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                             wing_semi_span_planform1=self.wing_semi_span_planform1,
                             wing_semi_span=self.corrected_semi_span,
                             wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                             wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                             front_spar_position=self.corrected_front_spar_position,
                             rear_spar_position=self.corrected_rear_spar_position,
                             rib_number=self.corrected_rib_number,

                             section_number=self.section_number,
                             points_number=self.points_number,
                             element_length=self.element_length
                             )

    @Attribute
    def material_choice(self):
        STEEL = DEFI_MATERIAU(ELAS=_F(E=300000000000.0, RHO=7850, NU=0.1666))
        ALUMINIUM = DEFI_MATERIAU(ELAS=_F(E=7e10, RHO=2700, NU=0.33))
        if self.material == 'Steel':
            return STEEL
        else:
            return ALUMINIUM

    @Attribute
    def fem_setup(self):
        return WingFEM(finalmesh=self.mesh,
                       avl=self.avl_analysis,
                       skin_writer=CodeAster_primitives(wing_airfoil_root=self.wing_airfoil_root,
                                                        wing_airfoil_middle=self.wing_airfoil_middle,
                                                        wing_airfoil_tip=self.wing_airfoil_tip,

                                                        wing_root_chord=self.corrected_root_chord,
                                                        wing_middle_chord=self.corrected_middle_chord,
                                                        wing_tip_chord=self.wing_tip_chord,

                                                        wing_thickness_factor_root=self.wing_thickness_factor_root,
                                                        wing_thickness_factor_middle=self.wing_thickness_factor_middle,
                                                        wing_thickness_factor_tip=self.wing_thickness_factor_tip,

                                                        wing_semi_span_planform1=self.wing_semi_span_planform1,
                                                        wing_semi_span=self.corrected_semi_span,
                                                        wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
                                                        wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,

                                                        front_spar_position=self.corrected_front_spar_position,
                                                        rear_spar_position=self.corrected_rear_spar_position,
                                                        rib_number=self.corrected_rib_number,

                                                        section_number=self.section_number,
                                                        points_number=self.points_number,

                                                        finalmesh=self.mesh
                                                        ))


    @Attribute
    def fem_writer(self):
        return Writer(instance=self.fem_setup, avl=self.avl_analysis, material=self.material_choice)

    @Attribute
    def run_fem_analysis(self):
        result = optimize_plate_thickness(
            target_deflection=self.target_deflection,
            wing_fem_instance=self.fem_setup,  # Pass your FEM instance
            writer_instance=self.fem_writer,  # Pass your Writer instance
            initial_thickness=self.corrected_initial_thickness,
            thickness_bounds=self.corrected_thickness_bounds
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
            wing_root_chord=self.corrected_root_chord,
            wing_middle_chord=self.corrected_middle_chord,
            wing_tip_chord=self.wing_tip_chord,
            wing_thickness_factor_root=self.wing_thickness_factor_root,
            wing_thickness_factor_middle=self.wing_thickness_factor_middle,
            wing_thickness_factor_tip=self.wing_thickness_factor_tip,
            wing_semi_span_planform1=self.wing_semi_span_planform1,
            wing_semi_span=self.corrected_semi_span,
            wing_sweep_leading_edge_planform1=self.wing_sweep_leading_edge_planform1,
            wing_sweep_leading_edge_planform2=self.wing_sweep_leading_edge_planform2,
            front_spar_thickness=float(self.run_fem_analysis['optimized_thickness']),
            front_spar_position=self.corrected_front_spar_position,
            rear_spar_thickness=float(self.run_fem_analysis['optimized_thickness']),
            rear_spar_position=self.corrected_rear_spar_position,

            rib_thickness=float(self.run_fem_analysis['optimized_thickness']),
            rib_number=self.corrected_rib_number,

            plate_thickness=float(self.run_fem_analysis['optimized_thickness']),

            stringer_thickness=self.stringer_thickness,
            stringer_number=self.stringer_number
        )

    # Step file creation here!!!
    @Part
    def step_writer(self):
        return STEPWriter(nodes=[self.wingbox.spars.solid_spar, self.wingbox.ribs, self.wingbox.plates])


if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
    )
    display(obj)
