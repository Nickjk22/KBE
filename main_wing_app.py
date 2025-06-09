from parapy.core import *
from parapy.gui import display
from parapy.geom import *
from wingbox import Wingbox
from AVL_analysis import WingAVLAnalysis
from FEM_analysis import WingFEM, Writer, ResultsReader
from torsionbox import TorsionBox
from visualisation_arrows import LiftArrowArray
import warnings
from FEM_analysis import optimize_plate_thickness
from wing import WingSurface

warnings.filterwarnings("ignore", category=UserWarning)  # Suppress AVL/FEM warnings


class IntegratedWingAnalysis(Base):
    # Wing Parameters
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

    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    rib_number = Input(15)

    stringer_thickness = Input(0.01)
    stringer_number = Input(10)

    # Results Storage
    # avl_results = Attribute()
    # fem_results = Attribute()
    # optimized_parameters = Attribute()

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

    # Stap 2: AVL-analyse (pas zichtbaar in GUI wanneer aangeklikt)
    @Attribute
    def avl_analysis(self):
        return WingAVLAnalysis(
            aircraft=self.wingbox.wing_surface,
            case_settings=[("alpha_5deg", {'alpha': 5.0})]
        )

    @Attribute
    def avl_lift_forces(self):
        lift_forces = []
        max_value = max(self.avl_analysis.lift_forces)
        for i in (self.avl_analysis.lift_forces):
            norm = i / max_value
            lift_forces.append(norm)
        return lift_forces

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

    @Part
    def lift_arrows(self):
        return LiftArrowArray(points_list=[pt.point for pt in TorsionBox(hidden=True).points],
                              lift_forces=self.avl_lift_forces)

    @Attribute
    def run_fem_analysis(self):
        result = optimize_plate_thickness(target_deflection=9)
        print(f"Optimized thickness: {result['optimized_thickness']} m"), print(
            f"Max deflection achieved: {result['max_deflection']} m"), print(
            f"Optimization success: {result['success']}")
        return result

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
            wing_twist=self.wing_twist,
            front_spar_thickness=self.run_fem_analysis['optimized_thickness'],
            front_spar_position=self.front_spar_position,
            rear_spar_thickness=self.run_fem_analysis['optimized_thickness'],
            rear_spar_position=Input(0.6),

            rib_thickness=self.run_fem_analysis['optimized_thickness'],
            rib_number=self.rib_number,

            stringer_thickness=self.stringer_thickness,
            stringer_number=self.stringer_number
        )

    # @Attribute
    # def design_report(self):
    #     """Shows after FEM completes"""
    #     report = (
    #         f"Optimized Design Parameters:\n"
    #         f"Front Spar Thickness: {self.optimized_parameters['front_spar_thickness']} m\n"
    #         f"Rear Spar Thickness: {self.optimized_parameters['rear_spar_thickness']} m\n"
    #         f"Rib Thickness: {self.optimized_parameters['rib_thickness']} m\n"
    #         f"Max Deflection: {self.fem_results.max_deflection} m"
    #     )
    #     # Make internals visible in wingbox
    #     self.wingbox.show_internals = True
    #     return report

if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
    )
    display(obj)
