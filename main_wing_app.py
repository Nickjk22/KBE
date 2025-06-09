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
    # avl_results = Attribute()
    # fem_results = Attribute()
    # optimized_parameters = Attribute()

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
            norm = i/max_value
            lift_forces.append(norm)
        return lift_forces

    @Attribute
    def torsionbox(self):
        return TorsionBox(hidden=True)

    @Attribute
    def points_list(self):
        return [pt.point for pt in self.torsionbox.points]

    # @Part
    # def lift_arrows(self):
    #     return Arrow(
    #         point=self.lift_arrow_data[child.index][0],
    #         direction=Vector(0, 0, 1) if self.lift_arrow_data[child.index][1] >= 0 else Vector(0, 0, -1),
    #         head_length=abs(self.lift_arrow_data[child.index][1]),
    #         base_length=abs(self.lift_arrow_data[child.index][1]),
    #         suppress=abs(self.lift_arrow_data[child.index][1]) < 1e-6
    #     ), Part.quantity(len(self.lift_arrow_data))

    @Part
    def lift_arrows(self):
        return LiftArrowArray(points_list=self.points_list, lift_forces=self.avl_lift_forces)

    # @Attribute
    # def fem_analysis(self):
    #     fem_instance = WingFEM(
    #         thickness=0.1,
    #     )
    #     writer = Writer(fem_instance)
    #
    #     # Write + Run
    #     import os
    #     base_path = os.path.dirname(__file__)
    #     output_dir = os.path.join(base_path, "output")
    #     os.makedirs(output_dir, exist_ok=True)
    #
    #     mesh_path = os.path.join(output_dir, "mesh2.aster")
    #     comm_path = os.path.join(output_dir, "wing_fem.comm")
    #     export_path = os.path.join(output_dir, "case2.export")
    #     results_path = os.path.join(output_dir, "results.txt")
    #     log_path = os.path.join(output_dir, "aster_log.txt")
    #
    #     writer.write_mesh(mesh_path)
    #     writer.write_comm(comm_path)
    #     from parapy.lib.code_aster import create_export_file, run_code_aster
    #     create_export_file(export_path, mesh_path, comm_path, results_path)
    #     run_code_aster(export_path, log_file=log_path)
    #
    #     return ResultsReader(file_path=results_path)
    #
    @Attribute
    def run_fem_analysis(self):
        result = optimize_plate_thickness(target_deflection=8.8)
        return print(f"Optimized thickness: {result['optimized_thickness']} m"), print(f"Max deflection achieved: {result['max_deflection']} m"), print(f"Optimization success: {result['success']}")


    # @Attribute
    # def avl_analysis(self):
    #     """Triggered when expanded in GUI tree"""
    #     if self.run_avl and not hasattr(self, 'avl_results'):
    #         self.avl_results = WingAVLAnalysis(
    #             aircraft=self.wingbox,
    #             case_settings=[("cruise", {'alpha': 5.0})]
    #         )
    #         return self.avl_results
    #     return "Run AVL by setting 'run_avl=True'"

    # @Attribute
    # def fem_analysis(self):
    #     """Triggered when expanded in GUI tree"""
    #     if self.run_fem and hasattr(self, 'avl_results'):
    #         if not hasattr(self, 'fem_results'):
    #             self.fem_results = WingFEM(
    #                 wing=self.wingbox,
    #                 avl=self.avl_results,
    #                 element_length=0.2,
    #                 thickness=0.1
    #             )
    #             # Update geometry with FEM results (example values)
    #             self.optimized_parameters = {
    #                 'front_spar_thickness': 0.07,
    #                 'rear_spar_thickness': 0.06,
    #                 'rib_thickness': 0.15,
    #                 'stringer_thickness': 0.02
    #             }
    #         return self.fem_results
    #     return "Run FEM by setting 'run_fem=True' (requires AVL first)"

    @Attribute
    def design_report(self):
        """Shows after FEM completes"""
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




if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
        run_avl=True,  # Set to True in GUI to run AVL
        run_fem=False  # Set to True in GUI to run FEM
    )
    display(obj)
