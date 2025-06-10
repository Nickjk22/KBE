from parapy.core import *
from parapy.gui import display
from parapy.geom import *
from wingbox import Wingbox
from AVL_analysis import WingAVLAnalysis
from visualisation_arrows import LiftArrowArray
import warnings
from FEM_analysis import optimize_plate_thickness
from wing import WingSurface
from meshing import FinalMesh
from FEM_analysis import WingFEM, Writer
from find_nodes import CodeAster_primitives

warnings.filterwarnings("ignore", category=UserWarning)  # Suppress AVL/FEM warnings


#Hier de airfoil coordinates interpoleren ipv torsionbox!!!!

class IntegratedWingAnalysis(Base):
    # Wing Parameters
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
    wing_twist = Input(0)

    front_spar_position = Input(0.2)
    rear_spar_position = Input(0.6)

    rib_number = Input(15)

    stringer_thickness = Input(0.01)
    stringer_number = Input(10)

    target_deflection = Input(5)
    initial_thickness = Input(0.3)
    thickness_bounds = Input([0.1, 0.5])

    # Results Storage
    # avl_results = Attribute()
    # fem_results = Attribute()
    # optimized_parameters = Attribute()
    points_number = Input(14)

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
                           wing_twist=self.wing_twist,
                           points_number=self.points_number
                           )

    # Stap 2: AVL-analyse (pas zichtbaar in GUI wanneer aangeklikt)
    @Attribute
    def avl_analysis(self):
        return WingAVLAnalysis(
            aircraft=self.wing_surface(label="wing"),
            case_settings=[("alpha_5deg", {'alpha': 5.0})],
            points_number=self.points_number,
            rho=1.2,
            Mach=0.7,
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
                                                        wing_twist=self.wing_twist,

                                                        front_spar_position=self.front_spar_position,
                                                        rear_spar_position=self.rear_spar_position,
                                                        rib_number=self.rib_number,

                                                        section_number=self.section_number,
                                                        segment_number=self.segment_number,
                                                        points_number=self.points_number,

                                                        finalmesh=self.mesh)),

        )

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

    @Part
    def lift_arrows(self):
        return LiftArrowArray(points_list=[pt.point for pt in self.wing_surface(hidden=True).points],
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
                         wing_twist=self.wing_twist,

                         front_spar_position=self.front_spar_position,
                         rear_spar_position=self.rear_spar_position,
                         rib_number=self.rib_number,

                         section_number=self.section_number,
                         segment_number=self.segment_number,
                         points_number=self.points_number)

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
                                                        wing_twist=self.wing_twist,

                                                        front_spar_position=self.front_spar_position,
                                                        rear_spar_position=self.rear_spar_position,
                                                        rib_number=self.rib_number,

                                                        section_number=self.section_number,
                                                        segment_number=self.segment_number,
                                                        points_number=self.points_number,

                                                        finalmesh=self.mesh))

    @Attribute
    def fem_writer(self):
        return Writer(





        )


#Hier inputs erbij van de 2 classes!
    @Attribute
    def run_fem_analysis(self):
        result = optimize_plate_thickness(target_deflection=self.target_deflection,
                                          initial_thickness=self.initial_thickness,
                                          thickness_bounds=self.thickness_bounds)
        print(f"Optimized thickness: {result['optimized_thickness']} m"), print(
            f"Max deflection achieved: {result['max_deflection']} m"), print(
            f"Optimization success: {result['success']}")
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
            wing_twist=self.wing_twist,
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


if __name__ == '__main__':
    obj = IntegratedWingAnalysis(
        label="Wing Analysis",
    )
    display(obj)
