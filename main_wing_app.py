from parapy.core import *
from parapy.gui import display
from parapy.geom import *
from wingbox import Wingbox
from AVL_analysis import WingAVLAnalysis
from FEM_analysis import WingFEM
import warnings

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

    # Analysis Control Flags
    run_avl = Input(False)
    run_fem = Input(False)

    # Results Storage
    avl_results = Attribute()
    fem_results = Attribute()
    optimized_parameters = Attribute()

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
            # Hide internals until FEM completes
            show_internals=hasattr(self, 'fem_results')
        )

    @Attribute
    def avl_analysis(self):
        """Triggered when expanded in GUI tree"""
        if self.run_avl and not hasattr(self, 'avl_results'):
            self.avl_results = WingAVLAnalysis(
                aircraft=self.wingbox,
                case_settings=[("cruise", {'alpha': 5.0})]
            )
            return self.avl_results
        return "Run AVL by setting 'run_avl=True'"

    @Attribute
    def fem_analysis(self):
        """Triggered when expanded in GUI tree"""
        if self.run_fem and hasattr(self, 'avl_results'):
            if not hasattr(self, 'fem_results'):
                self.fem_results = WingFEM(
                    wing=self.wingbox,
                    avl=self.avl_results,
                    element_length=0.2,
                    thickness=0.1
                )
                # Update geometry with FEM results (example values)
                self.optimized_parameters = {
                    'front_spar_thickness': 0.07,
                    'rear_spar_thickness': 0.06,
                    'rib_thickness': 0.15,
                    'stringer_thickness': 0.02
                }
            return self.fem_results
        return "Run FEM by setting 'run_fem=True' (requires AVL first)"

    @Attribute
    def design_report(self):
        """Shows after FEM completes"""
        if hasattr(self, 'fem_results'):
            report = (
                f"Optimized Design Parameters:\n"
                f"Front Spar Thickness: {self.optimized_parameters['front_spar_thickness']} m\n"
                f"Rear Spar Thickness: {self.optimized_parameters['rear_spar_thickness']} m\n"
                f"Rib Thickness: {self.optimized_parameters['rib_thickness']} m\n"
                f"Max Deflection: {self.fem_results.max_deflection} m"
            )
            # Make internals visible in wingbox
            self.wingbox.show_internals = True
            return report
        return "Complete FEM to see design report"


if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        run_avl=False,  # Set to True in GUI to run AVL
        run_fem=False  # Set to True in GUI to run FEM
    )
    display(obj)