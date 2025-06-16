from parapy.core import *
from parapy.exchange import STEPWriter
from parapy.gui import display
from wingbox import Wingbox
from AVL_analysis import WingAVLAnalysis
from visualisation_arrows import LiftArrowArray
import warnings
from FEM_analysis import optimize_plate_thickness
from wing import WingSurface
from meshing import MeshGenerator
from find_nodes import CodeAster_primitives
import numpy as np
import pandas as pd
import os
from scipy.interpolate import interp1d
from kbeutils import avl
import matplotlib.pyplot as plt
from parapy.lib.code_aster import (_F, DEFI_MATERIAU)
from parapy.core.validate import LessThanOrEqualTo, GreaterThan, GreaterThanOrEqualTo, Between, LessThan

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


# Usage
interpolate_airfoil('airfoil_data\whitcomb.dat', 'airfoil_data\whitcomb_interpolated.dat', factor=300)


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
    wing_airfoil_root = Input('airfoil_data\whitcomb_interpolated.dat')
    wing_airfoil_middle = Input("airfoil_data\whitcomb_interpolated.dat")
    wing_airfoil_tip = Input("airfoil_data\whitcomb_interpolated.dat")

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

    wing_sweep_leading_edge_planform1 = Input(float(pd.read_excel(excel_directory).iloc[11, 1]),
                                              validator=Between(-60, 60))
    wing_sweep_leading_edge_planform2 = Input(float(pd.read_excel(excel_directory).iloc[12, 1]),
                                              validator=Between(-60, 60))

    front_spar_position = Input(float(pd.read_excel(excel_directory).iloc[14, 1]), validator=Between(0, 1))
    rear_spar_position = Input(float(pd.read_excel(excel_directory).iloc[15, 1]), validator=Between(0, 1))

    @Attribute
    def corrected_front_spar_position(self):
        if self.front_spar_position + self.corrected_thickness_bounds2[
            1] / self.wing_tip_chord > self.rear_spar_position:
            msg = f"Front spar position ({self.front_spar_position}) plus its maximum thickness should be smaller than the rear spar position ({self.rear_spar_position}). Input will be ignored and front spar position will be placed to the front."
            warnings.warn(msg)
            if self.popup_gui_front_spar_position:
                generate_warning("Warning: Value changed", msg)
            return self.rear_spar_position - self.corrected_thickness_bounds2[1] / self.wing_tip_chord
        else:
            return self.front_spar_position

    @Attribute
    def corrected_rear_spar_position(self):
        if self.rear_spar_position + self.corrected_thickness_bounds2[1] / self.wing_tip_chord > 1:
            msg = f"Rear spar position ({self.rear_spar_position}) plus its maximum thickness should be smaller than the chord length. Input will be ignored and rear spar position will be placed to the front."
            warnings.warn(msg)
            if self.popup_gui_rear_spar_position:
                generate_warning("Warning: Value changed", msg)
            return 1 - self.corrected_thickness_bounds2[1] / self.wing_tip_chord
        else:
            return self.rear_spar_position

    stringer_thickness = Input(float(pd.read_excel(excel_directory).iloc[19, 1]), validator=GreaterThan(0))
    stringer_number = Input(int(pd.read_excel(excel_directory).iloc[20, 1]), validator=GreaterThan(0))

    target_deflection = Input(float(pd.read_excel(excel_directory).iloc[23, 5]), validator=GreaterThan(0))
    thickness_bounds = Input(
        [float(pd.read_excel(excel_directory).iloc[25, 5]), float(pd.read_excel(excel_directory).iloc[26, 5])])
    initial_thickness = Input(float(pd.read_excel(excel_directory).iloc[24, 5]), validator=GreaterThan(0))

    @Attribute
    def airfoil_thickness(self):
        # Pad naar het bestand
        path = self.wing_airfoil_root if isinstance(self.wing_airfoil_root, str) else self.wing_airfoil_root.path

        # Read data from file
        with open(path, 'r') as f:
            coords = [tuple(map(float, line.strip().split())) for line in f if line.strip()]

        # Extract X and Y separately
        x_vals = [pt[0] for pt in coords]

        # Search the x-coördinate that is the closest to front_spar_position
        x_front = min(x_vals, key=lambda x: abs(x - self.front_spar_position))
        x_rear = min(x_vals, key=lambda x: abs(x - self.rear_spar_position))

        # Extract the y-values on both sides of the x value
        front_ys = [y for x, y in coords if abs(x - x_front) < 1e-5]
        rear_ys = [y for x, y in coords if abs(x - x_rear) < 1e-5]

        # Check if the x-values have an upper and lower side
        if len(front_ys) < 2 or len(rear_ys) < 2:
            raise ValueError("Niet genoeg y-punten gevonden bij sparlocaties.")

        # sort y values to sort by upper and lower
        front_ys.sort()
        rear_ys.sort()

        # Take the difference between upper and lower
        front_thickness = abs(front_ys[-1] - front_ys[0])
        rear_thickness = abs(rear_ys[-1] - rear_ys[0])

        # Choose the smallest value for the thickness
        return min(front_thickness, rear_thickness)

    @Attribute
    def corrected_thickness_bounds1(self):
        if self.thickness_bounds[1] > self.airfoil_thickness * self.wing_tip_chord * self.wing_thickness_factor_tip:
            msg = (f"Maximum plate thickness ({self.thickness_bounds[1]}) should be smaller than the thinnest part of "
                   f"the airfoil between the spars. Input will be ignored and upper bound will be set equal to the "
                   f"thinnest part")
            warnings.warn(msg)
            if self.popup_gui_thickness_bounds:
                generate_warning("Warning: Value changed", msg)
            return [self.thickness_bounds[0],
                    (self.airfoil_thickness * self.wing_tip_chord * self.wing_thickness_factor_tip)]
        else:
            return [self.thickness_bounds[0], self.thickness_bounds[1]]

    @Attribute
    def corrected_thickness_bounds2(self):
        if self.corrected_thickness_bounds1[0] > self.corrected_thickness_bounds1[1]:
            msg = (f"Lower thickness bound ({self.c1[0]}) should be smaller than the upper thickness "
                   f"bound ({self.corrected_thickness_bounds1[1]}). Input will be ignored and lower bound will be set equal to "
                   f"upper bound.")
            warnings.warn(msg)
            if self.popup_gui_thickness_bounds:
                generate_warning("Warning: Value changed", msg)
            return [self.corrected_thickness_bounds1[1], self.corrected_thickness_bounds1[1]]
        else:
            return [self.corrected_thickness_bounds1[0], self.corrected_thickness_bounds1[1]]

    @Attribute
    def corrected_initial_thickness(self):
        if not self.corrected_thickness_bounds2[0] < self.initial_thickness < self.corrected_thickness_bounds2[1]:
            msg = f"Initial thickness ({self.initial_thickness}) should be in between lower and upper bound ({self.corrected_thickness_bounds2}). Input will be ignored and set to average of bounds."
            warnings.warn(msg)
            if self.popup_gui_initial_thickness:
                generate_warning("Warning: Value changed", msg)
            return (self.corrected_thickness_bounds2[0] + self.corrected_thickness_bounds2[1]) / 2
        else:
            return self.initial_thickness

    rib_number = Input(int(pd.read_excel(excel_directory).iloc[17, 1]), validator=GreaterThan(0))

    @Attribute
    def corrected_rib_number(self):
        if self.rib_number * self.corrected_thickness_bounds2[1] > self.corrected_semi_span:
            msg = f"Number of ribs ({self.rib_number}) times maximum rib thickness (thickness upper bound) ({self.corrected_thickness_bounds2[1]}) exceeds the wing semi-span. Input will be ignored and the rib number will be set equal to wing semi-span divided by the maximum thickness."
            warnings.warn(msg)
            if self.popup_gui_rib_number:
                generate_warning("Warning: Value changed", msg)
            return self.wing_semi_span / self.corrected_thickness_bounds2[1]
        else:
            return self.rib_number

    points_number = Input(int(pd.read_excel(excel_directory).iloc[22, 1]), validator=GreaterThan(0))
    section_number = Input(int(pd.read_excel(excel_directory).iloc[23, 1]), validator=GreaterThan(0))

    Mach = Input(float(pd.read_excel(excel_directory).iloc[17, 5]), validator=Between(0, 1))
    rho = Input(float(pd.read_excel(excel_directory).iloc[20, 5]), validator=GreaterThan(0))

    element_length = Input(float(pd.read_excel(excel_directory).iloc[27, 5]), validator=Between(0, 1))
    is_mirrored = Input(True)

    material = Input(str(pd.read_excel(excel_directory).iloc[28, 5]))

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

    # Read the settings based on case number (1 = fixed CL, 2 = fixed angle of attack)
    case_setting = Input(float(pd.read_excel(excel_directory).iloc[18, 5]))
    new_parameter = Input(float(pd.read_excel(excel_directory).iloc[19, 5]))

    # Write the cases based on the input
    @Attribute
    def case_inputs(self):
        if self.case_setting == 2:
            case1 = ('fixed_aoa', {'alpha': self.new_parameter})
            return case1
        else:
            case2 = ('fixed_cl', {'alpha': avl.Parameter(name='alpha', value=self.new_parameter, setting='CL')})
            return case2

    # Stap 2: AVL-analysis
    @Attribute
    def avl_analysis(self):
        return WingAVLAnalysis(
            aircraft=self.wing_surface
            # (label="wing")
            ,
            case_settings=[self.case_inputs],
            is_mirrored=self.is_mirrored,
            points_number=(2 * self.points_number) if self.is_mirrored else self.points_number,
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

                                             finalmesh=self.mesh,
                                             element_length=self.element_length))

    @Attribute
    def avl_and_optimisation_results(self):
        new_file_name = "avl_and_optimisation_results.txt"
        directory_path = excel_directory.rsplit('\\', 1)[0]
        file_path = directory_path + '\\' + 'Project\\KBE\\output' + '\\' + new_file_name
        result = self.run_fem_analysis

        # Clear the output file and rewrite content, including the inputs used for AVL, geometry and optimisation
        with open(file_path, 'w') as f:
            f.write("Output file for the AVL and FEM analysis")
            f.write('\n')
            f.write('\n')
            f.write('AVL inputs:')
            f.write('\n')
            f.write(f"Fixed aoa or Cl: {self.new_parameter} m")
            f.write('\n')
            f.write(f"Mach number: {self.Mach} m")
            f.write('\n')
            f.write(f"Air density: {self.rho} m")
            f.write('\n')
            f.write('\n')
            f.write('Inputs used in this simulation:')
            f.write('\n')
            f.write(f"Wing root chord: {self.wing_root_chord} m")
            f.write('\n')
            f.write(f"Wing middle chord: {self.wing_middle_chord} m")
            f.write('\n')
            f.write(f"Wing tip chord: {self.wing_tip_chord} m")
            f.write('\n')
            f.write(f"Wing thickness factor root: {self.wing_thickness_factor_root} m")
            f.write('\n')
            f.write(f"Wing thickness factor middle: {self.wing_thickness_factor_middle} m")
            f.write('\n')
            f.write(f"Wing thickness factor tip: {self.wing_thickness_factor_tip} m")
            f.write('\n')
            f.write(f"Wing semi span: {self.wing_semi_span} m")
            f.write('\n')
            f.write(f"Wing span planform 1: {self.wing_semi_span_planform1} m")
            f.write('\n')
            f.write(f"Wing sweep of leading edge of planform 1: {self.wing_sweep_leading_edge_planform1} deg")
            f.write('\n')
            f.write(f"Wing sweep of leading edge of planform 2: {self.wing_sweep_leading_edge_planform2} deg")
            f.write('\n')
            f.write(f"Front spar position: {self.front_spar_position} m")
            f.write('\n')
            f.write(f"Rear spar position: {self.rear_spar_position} m")
            f.write('\n')
            f.write(f"Stringer number: {self.stringer_number} m")
            f.write('\n')
            f.write(f"Rib number: {self.rib_number} m")
            f.write('\n')
            f.write(f"Stringer thickness: {self.stringer_thickness} m")
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write(f"Target deflection: {self.target_deflection} m")
            f.write('\n')
            f.write(f"Thickness bounds: [{self.thickness_bounds[0]}, {self.thickness_bounds[1]}] m")
            f.write('\n')
            f.write(f"Initial thickness: {self.initial_thickness} m")
            f.write('\n')
            f.write(f"Material: {self.material} m")
            f.write('\n')
            f.write('\n')
            f.write('\n')
            f.write(f"Optimized thickness: {result['optimized_thickness']} m")
            f.write('\n')
            f.write(f"Max deflection achieved: {result['max_deflection']} m")
            f.write('\n')
            f.write(f"Optimization success: {result['success']}")
            f.write('\n')
            f.write('\n')

        with open(file_path, 'a') as f:
            f.write("Lift Forces per point:\n")
            for i, value in enumerate(self.avl_lift_forces):
                f.write(f"Lift at point {i}: {value:.6f} N\n")
            f.write('\n')

        # Print per case the totals and the stability derivatives
        for case_name, result in self.avl_analysis.results.items():
            with open(file_path, 'a') as f:
                f.write("Case name = ")
                f.write(str(result['Name']))
                f.write('\n')
                f.write('\n')

                # Totals
                Totals = ['CYtot', 'Cmtot', 'CZtot', 'Cntot', 'CLtot', 'CDtot', 'CDind', 'CLff']
                f.write('Totals:')
                f.write('\n')
                for i in Totals:
                    f.write(i)
                    f.write(';   ')
                    f.write(str(result['Totals'][i]))
                    f.write('\n')
                f.write('\n')

                # Stability derivatives
                Stability_der = ['CLa', 'CLb', 'CYa', 'CYb', 'Cla', 'Clb', 'Cma', 'Cmb', 'Cna', 'Cnb', 'CLp', 'CLq',
                                 'CLr', 'CYp', 'CYq', 'CYr', 'Clp', 'Clq', 'Clr', 'Cmp', 'Cmq', 'Cmr', 'Cnp', 'Cnq',
                                 'Cnr', 'Xnp']
                f.write('Stability derivatives:')
                f.write('\n')
                for i in Stability_der:
                    f.write(i)
                    f.write(';   ')
                    f.write(str(result['StabilityDerivatives'][i]))
                    f.write('\n')

        # return self.results.items()
        return "results generated in text file"

    @Attribute
    def plot_lift_distribution(self):
        lift_values = self.avl_lift_forces
        num_points = len(lift_values)

        # Generate evenly spaced spanwise positions from 0 to wing_semi_span
        spanwise_positions = np.linspace(0, self.wing_semi_span, num_points)

        # Plot
        plt.figure(figsize=(10, 5))
        plt.plot(spanwise_positions, lift_values, label="Lift Distribution", color="blue")
        plt.xlabel("Spanwise Position [m]")
        plt.ylabel("Lift Force [N]")
        plt.title("Lift Distribution over Semi-Span")
        plt.grid(True)
        plt.legend()
        plt.xlim(0, self.wing_semi_span)

        # Save the figure as .eps in the output folder
        output_dir = os.path.join(os.path.dirname(__file__), "output")
        os.makedirs(output_dir, exist_ok=True)

        file_path = os.path.join(output_dir, "lift_distribution.eps")
        plt.savefig(file_path, format='eps')
        plt.close()

        return file_path

    @Attribute
    def avl_lift_forces_normalized(self):
        lift_forces = []
        max_value = max(self.avl_analysis.lift_forces)
        for i in self.avl_analysis.lift_forces:
            norm = i / max_value
            lift_forces.append(norm)
        return lift_forces

    @Attribute
    def avl_lift_forces(self):
        return self.avl_analysis.lift_forces

    @Attribute
    def points_list(self):
        return self.wing_surface.points  # Points list based on points_number

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
        STEEL = DEFI_MATERIAU(ELAS=_F(E=3e11, RHO=7850, NU=0.1666))  # structural steel
        ALUMINIUM = DEFI_MATERIAU(ELAS=_F(E=7e10, RHO=2700, NU=0.33))  # generic aluminum
        AL7075 = DEFI_MATERIAU(ELAS=_F(E=7.1e10, RHO=2810, NU=0.33))  # Aluminum 7075-T6
        TI6AL4V = DEFI_MATERIAU(ELAS=_F(E=1.1e11, RHO=4430, NU=0.34))  # Titanium alloy
        CFRP = DEFI_MATERIAU(ELAS=_F(E=1.4e11, RHO=1600, NU=0.28))  # Carbon Fiber Reinforced Polymer (average)
        GFRP = DEFI_MATERIAU(ELAS=_F(E=3.5e10, RHO=1900, NU=0.27))  # Glass Fiber Reinforced Polymer

        if self.material == 'Steel':
            return STEEL
        elif self.material == 'Aluminium':
            return ALUMINIUM
        elif self.material == 'AL7075':
            return AL7075
        elif self.material == 'Titanium':
            return TI6AL4V
        elif self.material == 'CFRP':
            return CFRP
        elif self.material == 'GFRP':
            return GFRP
        else:
            raise ValueError(f"Unknown material: {self.material}")

    @Attribute
    def run_fem_analysis(self):  # Runs optimisation
        result = optimize_plate_thickness(
            target_deflection=self.target_deflection,
            initial_thickness=self.corrected_initial_thickness,
            thickness_bounds=self.corrected_thickness_bounds2,
            check_element=0,
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
            element_length=self.element_length,
            avl_analysis=self.avl_analysis,
            material_choice=self.material_choice,
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

                                             finalmesh=self.mesh,
                                             element_length=self.element_length
                                             )
        )
        print(f"Optimized thickness: {result['optimized_thickness']} m")
        print(f"Max deflection achieved: {result['max_deflection']} m")
        print(f"Optimization success: {result['success']}")
        return result

    # Showcase optimized thicknesses, and generate stepfile
    @Part()
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

    @Part
    def step_writer(self):
        """
        Generates a STEPWriter object for exporting wingbox components.
        """
        return STEPWriter(nodes=[self.wingbox.spars.front_spar,
                                 self.wingbox.spars.rear_spar,
                                 self.wingbox.plates.upper_plate,
                                 self.wingbox.plates.lower_plate] +
                                [rib.lofted_solid for rib in self.wingbox.ribs])

    @Attribute
    def wing_mass(self):
        total_volume = (self.wingbox.spars.volume +
                        self.wingbox.plates.volume +
                        sum(rib.volume for rib in self.wingbox.ribs))
        density = self.material_choice.ELAS.RHO
        return total_volume * density


if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
    )
    display(obj)
