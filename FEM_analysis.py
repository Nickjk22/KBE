# -*- coding: utf-8 -*-
#
# Copyright (C) 2016-2023 ParaPy Holding B.V.
#
# You may use the contents of this file in your application code.
#
# THIS CODE AND INFORMATION ARE PROVIDED "AS IS" WITHOUT WARRANTY OF ANY
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR
# PURPOSE.

"""Demonstrates how to set up a 2D FEM analysis. Constraining is done by
clamping half of the nodes on the sides. Loading is applied via point loads.
"""

from collections import namedtuple
from typing import Dict, List, NamedTuple
from scipy.optimize import minimize
import os

from parapy.core import Attribute, Base, Input, Part, child
from glueing import GeneralFuse
from parapy.mesh import EdgeGroup, FaceGroup
from parapy.mesh.salome import Mesh, Tri
from parapy.geom import Compound, Face
from parapy.mesh.salome import MeshNode
from parapy.mesh.salome.controls import FixedLength
from parapy.lib.code_aster import (_F, AFFE_CARA_ELEM, AFFE_CHAR_MECA,
                                   AFFE_MATERIAU, AFFE_MODELE, DEFI_MATERIAU,
                                   IMPR_RESU, LIRE_MAILLAGE, MECA_STATIQUE,
                                   Command, CommandWriter, MeshGroup,
                                   MeshWriter, ResultsReaderBase,
                                   create_export_file, run_code_aster)
from meshing import MeshGenerator
from AVL_analysis import WingAVLAnalysis
from find_nodes import CodeAster_primitives
from wing import WingSurface

# FACE = "face_group"
CONSTRAINED_EDGE1 = "constrained_edge1_group"
CONSTRAINED_EDGE2 = "constrained_edge2_group"
CONSTRAINED_EDGE3 = "constrained_edge3_group"
CONSTRAINED_EDGE4 = "constrained_edge4_group"


# LOADED_EDGE1 = "loaded_edge1_group"
# LOADED_EDGE2 = "loaded_edge2_group"


class WingFEM(Base):
    thickness: float = Input(0.1)
    length: float = Input(1)
    width: float = Input(2)
    element_length: float = Input(0.1)

    @Input
    def finalmesh(self):
        return MeshGenerator(element_length=self.element_length)

    @Attribute
    def shape_to_mesh(self) -> GeneralFuse:
        return self.finalmesh.shape_to_mesh

    @Part
    def contrained_edge_groups(self) -> EdgeGroup:
        return EdgeGroup(quantify=4,
                         shape=[self.shape_to_mesh.edges[3],
                                self.shape_to_mesh.edges[6],
                                self.shape_to_mesh.edges[9],
                                self.shape_to_mesh.edges[10]][child.index],
                         label=[CONSTRAINED_EDGE1,
                                CONSTRAINED_EDGE2,
                                CONSTRAINED_EDGE3,
                                CONSTRAINED_EDGE4][child.index])

    # @Part
    # def loaded_edge_group(self) -> EdgeGroup:
    #     return EdgeGroup(quantify=2,
    #                      shape=[self.shape_to_mesh.edges[156], self.shape_to_mesh.edges[157]][child.index],
    #                      label=[LOADED_EDGE1, LOADED_EDGE2][child.index])

    @Attribute
    def FACE(self):
        return [f"face{i + 1}_group" for i in range(len(self.shape_to_mesh.faces))]

    @Part
    def face_group(self) -> FaceGroup:
        return FaceGroup(quantify=len(self.shape_to_mesh.faces),
                         shape=[self.shape_to_mesh.faces[child.index]],
                         label=self.FACE[child.index])

    @Attribute
    def mesh(self) -> Mesh:
        return self.finalmesh.mesh

    # @Input
    # def nodes(self) -> List[int]:
    #     """List of mesh‐node IDs to apply lift forces on."""
    #     return []
    #
    # @Input
    # def liftforces(self) -> List[float]:
    #     """Same length as `nodes`: force magnitude (FZ) per node."""
    #     return []

    @Input
    def avl(self):
        return WingAVLAnalysis(aircraft=WingSurface(label="wing"), case_settings=[
            ("alpha_5deg", {'alpha': 5.0}),
        ])

    @Input
    def skin_writer(self):
        return CodeAster_primitives()

    @Attribute
    def liftforces(self) -> List[float]:
        return self.avl.lift_forces

    @Attribute
    def nodes(self) -> List[MeshNode]:
        return self.skin_writer.load_primitives


class Writer:
    """Writes code_aster command and mesh files.

    # >>> instance = WingFEM()
    # >>> obj = Writer(instance)
    # >>> obj.write_comm("C:/Users/raane/Documents/Uni/Master/KBE/Year2/Tutorials/plate_CodeAster/output/output.comm")  # doctest: +ELLIPSIS
    # # >>> obj.write_comm("C:/Users/nick2/PycharmProjects//KBE/GitHub/output/output.comm")  # doctest: +ELLIPSIS
    #
    # # Written: ...
    # >>> obj.write_mesh("C:/Users/raane/Documents/Uni/Master/KBE/Year2/Tutorials/plate_CodeAster/output/mesh2.aster")  # doctest: +ELLIPSIS
    # # >>> obj.write_mesh("C:/Users/nick2/PycharmProjects//KBE/GitHub/output/mesh2.aster")  # doctest: +ELLIPSIS

    Written: ...
    """

    def __init__(self, instance: WingFEM = None, avl: WingAVLAnalysis = None, material: DEFI_MATERIAU = None) -> None:
        self._instance: WingFEM = instance or self._default_instance()
        self.avl: WingAVLAnalysis = avl or self._default_avl()
        self.material: DEFI_MATERIAU = material or self._default_material()

        self.mesh_settings_command: LIRE_MAILLAGE = None
        self.mesh_groups: List[MeshGroup] = None
        self.model_command: AFFE_MODELE = None
        self.shell_properties_command: AFFE_CARA_ELEM = None
        self.material_zone_command: AFFE_MATERIAU = None
        self.constraint_commands: List[AFFE_CHAR_MECA] = None
        self.load_commands: List[AFFE_CHAR_MECA] = None
        self.solver_settings_command: MECA_STATIQUE = None
        self.result_settings_command: IMPR_RESU = None

        self._generate_mesh_groups()
        self._generate_mesh_settings_command()
        self._generate_model_command()
        self._generate_shell_properties_command()
        self._generate_material_zone_command()
        self._generate_constraint_commands()
        self._generate_load_commands()
        self._generate_solver_settings_command()
        self._generate_results_settings_command()

        # Commands that are part of other commands do not need to be
        # collected. These nested commands are unpacked recursively by the
        # Code Aster writer. Examples of nested commands are:
        #   MAILLAGE in AFFE_MODELE
        #   AFFE in AFFE_MODELE, etc
        # As such, we only need to pass the root command which recursively
        # contains all other commands of the FEM model.
        self._commands: List[Command] = [self.result_settings_command]

    # @Part
    # def skin_writer(self):
    #     return CodeAster_primitives()

    @staticmethod
    def _default_instance() -> WingFEM:
        """Fallback WingFEM setup for standalone execution."""
        # Example placeholder — replace with actual logic or test setup
        return WingFEM()

    @staticmethod
    def _default_avl() -> WingAVLAnalysis:
        """Fallback AVL setup for standalone execution."""
        # You’ll likely need to provide a minimal working AVL configuration
        wing = WingSurface(label="wing")
        return WingAVLAnalysis(aircraft=wing, case_settings=[("alpha_5deg", {'alpha': 5.0})])

    @staticmethod
    def _default_material() -> DEFI_MATERIAU:
        ALUMINIUM = DEFI_MATERIAU(ELAS=_F(E=7e10, RHO=2700, NU=0.33))
        return ALUMINIUM

    @Attribute
    def liftforces(self):
        return self.avl.lift_forces

    def write_comm(self, pathname: str) -> None:
        writer = CommandWriter(self._commands)
        return writer.write(pathname)

    def write_mesh(self, pathname: str) -> None:
        writer = MeshWriter(grid=self._instance.mesh.grid,
                            groups=self.mesh_groups)
        return writer.write(pathname)

    def _generate_mesh_groups(self) -> None:
        instance = self._instance
        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE1,
                                                        shape=self._instance.shape_to_mesh.edges[3]).nodes
        ids = [elm.mesh_id for elm in elements]
        group1 = MeshGroup(label=CONSTRAINED_EDGE1,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE2,
                                                        shape=self._instance.shape_to_mesh.edges[6]).nodes
        ids = [elm.mesh_id for elm in elements]
        group2 = MeshGroup(label=CONSTRAINED_EDGE2,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE3,
                                                        shape=self._instance.shape_to_mesh.edges[9]).nodes
        ids = [elm.mesh_id for elm in elements]
        group3 = MeshGroup(label=CONSTRAINED_EDGE3,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE4,
                                                        shape=self._instance.shape_to_mesh.edges[10]).nodes
        ids = [elm.mesh_id for elm in elements]
        group4 = MeshGroup(label=CONSTRAINED_EDGE4,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        # elements = instance.mesh.get_subgrid_on_the_fly(label=LOADED_EDGE1,
        #                                                 shape=self._instance.shape_to_mesh.edges[157]).nodes
        # ids = [elm.mesh_id for elm in elements]
        # group5 = MeshGroup(label=LOADED_EDGE1,
        #                    header="GROUP_NO",
        #                    element_ids=ids,
        #                    element_type="node")

        # elements = instance.mesh.get_subgrid_on_the_fly(label=LOADED_EDGE2,
        #                                                 shape=self._instance.shape_to_mesh.edges[156]).nodes
        # ids = [elm.mesh_id for elm in elements]
        # group6 = MeshGroup(label=LOADED_EDGE2,
        #                    header="GROUP_NO",
        #                    element_ids=ids,
        #                    element_type="node")

        self.mesh_groups = [group1, group2, group3, group4,
                            # group5, group6
                            ]

        self.FACE = [f"face{i + 1}_group" for i in range(len(self._instance.shape_to_mesh.faces))]

        # Create and append MeshGroups for each face starting from group6
        for i, face_label in enumerate(self.FACE):
            shape = self._instance.shape_to_mesh
            element_ids = []

            if isinstance(shape, Compound):
                # Get all subfaces directly from the compound
                faces = shape.faces
            elif isinstance(shape, Face):
                faces = [shape]
            else:
                # fallback: try to access `.faces` if possible
                faces = getattr(shape, 'faces', [])

            for face in faces:
                try:
                    subgrid = instance.mesh.get_subgrid_on_the_fly(
                        label=face_label,
                        shape=face
                    )
                    element_ids.extend(f.mesh_id for f in subgrid.faces)
                except Exception as e:
                    print(f"Warning: Failed to process face for label {face_label}: {e}")

            group = MeshGroup(
                label=face_label,
                header="GROUP_MA",
                element_ids=element_ids,
                element_type="face"
            )
            self.mesh_groups.append(group)

        for i, mesh_node in enumerate(self._instance.nodes):
            grp_label = f"load_node{i + 1}_group"
            node_group = MeshGroup(
                label=grp_label,
                header="GROUP_NO",
                element_ids=[mesh_node.mesh_id],  # directly use the mesh_id
                element_type="node"
            )
            self.mesh_groups.append(node_group)

    def _generate_mesh_settings_command(self) -> None:
        self.mesh_settings_command = LIRE_MAILLAGE(UNITE=20,
                                                   FORMAT="ASTER")

    def _generate_model_command(self) -> None:
        shell_model = _F(GROUP_MA=tuple(self.FACE),
                         MODELISATION=("DKT",),
                         PHENOMENE="MECANIQUE")
        self.model_command = AFFE_MODELE(AFFE=(shell_model,),
                                         MAILLAGE=self.mesh_settings_command)

    def _generate_shell_properties_command(self) -> None:
        # Reference vector needed by CodeAster is based on normal and first
        # tangent at first node of first face of subgrid.
        # see CodeAster docs U4.42.01:8.3.4 for more info.
        face = self._instance.mesh.grid.faces[0]
        tangent = face.tangent1
        normal = face.plane_normal
        vec = tangent + normal
        self.shell_properties_command = AFFE_CARA_ELEM(
            MODELE=self.model_command,
            COQUE=[_F(
                EPAIS=self._instance.thickness,
                GROUP_MA=tuple(self.FACE),
                VECTEUR=(vec.x, vec.y, vec.z))])

    def _generate_material_zone_command(self) -> None:
        self.material_zone_command = AFFE_MATERIAU(AFFE=(_F(GROUP_MA=tuple(self.FACE),
                                                            MATER=(self.material,)),),
                                                   MAILLAGE=self.mesh_settings_command,
                                                   MODELE=self.model_command)

    def _generate_constraint_commands(self) -> None:
        self.constraint_commands = [AFFE_CHAR_MECA(
            DDL_IMPO=_F(
                GROUP_NO=(CONSTRAINED_EDGE1, CONSTRAINED_EDGE2, CONSTRAINED_EDGE3, CONSTRAINED_EDGE4),
                LIAISON="ENCASTRE"),
            MODELE=self.model_command)]

    # def _generate_load_commands(self) -> None:
    #     self.load_commands = [AFFE_CHAR_MECA(FORCE_NODALE=_F(FZ=40,
    #                                                          GROUP_NO=(LOADED_EDGE1,)),
    #                                          MODELE=self.model_command), AFFE_CHAR_MECA(FORCE_NODALE=_F(FZ=1,
    #                                                                                                     GROUP_NO=(
    #                                                                                                         LOADED_EDGE2,)),
    #                                                                                     MODELE=self.model_command)]

    def _generate_load_commands(self) -> None:
        # print(">>> DEBUG: nodes.mesh_id  =", [n.mesh_id for n in self._instance.nodes])
        # print(">>> DEBUG: liftforces      =", self._instance.liftforces)
        self.load_commands = []
        for i, force in enumerate(self._instance.liftforces):
            grp_label = f"load_node{i + 1}_group"
            cmd = AFFE_CHAR_MECA(
                FORCE_NODALE=_F(
                    GROUP_NO=(grp_label,),
                    FZ=force,
                ),
                MODELE=self.model_command
            )
            self.load_commands.append(cmd)

    def _generate_solver_settings_command(self) -> None:
        loads_constraints = [_F(CHARGE=obj) for obj in self.load_commands + self.constraint_commands]
        self.solver_settings_command = MECA_STATIQUE(SOLVEUR=_F(RESI_RELA=1e-06),
                                                     CARA_ELEM=self.shell_properties_command,
                                                     CHAM_MATER=self.material_zone_command,
                                                     EXCIT=loads_constraints,
                                                     MODELE=self.model_command)

    def _generate_results_settings_command(self):
        self.result_settings_command = IMPR_RESU(FORMAT="RESULTAT",
                                                 RESU=(_F(NOM_CHAM=("DEPL",),
                                                          RESULTAT=self.solver_settings_command),),
                                                 UNITE=8)


class ResultsReader(ResultsReaderBase):
    @Attribute
    def nodes_displacements(self) -> Dict[str, NamedTuple]:
        results_ = self.results["DEPL"]
        tmplt = namedtuple("nodes_displacements", ["TX", "TY", "TZ", "RX", "RY", "RZ"])
        dct = {}
        for line in results_[1:]:
            line_ = line.split()
            dct[line_[0]] = tmplt(*line_[1:])

        return dct

    @Attribute
    def max_deflection(self) -> float:
        deflections = [float(value.TZ) for value in self.nodes_displacements.values()]
        return max(deflections)


# if __name__ == "__main__":
#     instance = WingFEM(finalmesh=MeshGenerator(
#                        check_element=0,
#                        wing_airfoil_root="whitcomb_interpolated.dat",
#                        wing_airfoil_middle="whitcomb_interpolated.dat",
#                        wing_airfoil_tip="whitcomb_interpolated.dat",
#
#                        wing_root_chord=6,
#                        wing_middle_chord=4,
#                        wing_tip_chord=1.5,
#
#                        wing_thickness_factor_root=1,
#                        wing_thickness_factor_middle=1,
#                        wing_thickness_factor_tip=1,
#
#                        wing_semi_span_planform1=5,
#                        wing_semi_span=16,
#                        wing_sweep_leading_edge_planform1=20,
#                        wing_sweep_leading_edge_planform2=20,
#
#                        front_spar_position=0.2,
#                        rear_spar_position=0.6,
#                        rib_number=12,
#
#                        section_number=14,
#                        points_number=14,
#                        element_length=0.1
#                        ), element_length=0.1
#     )
#
#     current_path = os.path.dirname(__file__)
#     output_dir = os.path.join(current_path, "output")
#
#     writer = Writer(instance)
#     mesh_path = os.path.join(output_dir, "mesh2.aster")
#     writer.write_mesh(mesh_path)
#
#     comm_path = os.path.join(output_dir, "wing_fem.comm")
#     writer.write_comm(comm_path)
#
#     export_path = os.path.join(output_dir, "case2.export")
#     results_path = os.path.join(output_dir, "results.txt")
#     create_export_file(
#         target_path=export_path,
#         mesh_path=mesh_path,
#         comm_path=comm_path,
#         results_path=results_path,
#     )
#
#     log_file = os.path.join(output_dir, "aster_log.txt")
#     run_code_aster(export_path, log_file=log_file)
#
#     results_reader = ResultsReader(file_path=results_path)
#     print(f"maximum deflection from FEM: {float(results_reader.max_deflection)}")
#


# def optimize_plate_thickness(target_deflection: float,
#                              wing_fem_instance,  # New argument for WingFEM instance
#                              writer_instance,  # New argument for Writer instance
#                              initial_thickness=0.1,
#                              thickness_bounds=(0.02, 0.5)):
#     """
#     Optimize plate thickness to minimize thickness while keeping max deflection <= target_deflection.
#
#     Args:
#         target_deflection (float): Maximum allowed deflection (e.0.01 meters).
#         wing_fem_instance: An instance of the WingFEM class.
#         writer_instance: An instance of the Writer class.
#         initial_thickness (float): Initial guess for the plate thickness.
#         thickness_bounds (tuple): (min_thickness, max_thickness) allowed bounds.
#
#     Returns:
#         dict: Contains optimized thickness, max deflection, success status, and optimizer info.
#     """
#
#     def objective(thickness):
#         """Objective: minimize thickness."""
#         return thickness[0]  # minimize thickness itself
#
#     def constraint(thickness):
#         """Constraint: deflection must be <= target_deflection."""
#         # --- Update WingFEM instance with new thickness ---
#         # wing_fem_instance.thickness = thickness[0]
#
#         print("Thickness =", thickness[0])
#
#         wing_fem_instance = WingFEM(finalmesh=MeshGenerator(
#             check_element=0,
#             wing_airfoil_root="whitcomb_interpolated.dat",
#             wing_airfoil_middle="whitcomb_interpolated.dat",
#             wing_airfoil_tip="whitcomb_interpolated.dat",
#
#             wing_root_chord=6,
#             wing_middle_chord=4,
#             wing_tip_chord=1.5,
#
#             wing_thickness_factor_root=1,
#             wing_thickness_factor_middle=1,
#             wing_thickness_factor_tip=1,
#
#             wing_semi_span_planform1=5,
#             wing_semi_span=16,
#             wing_sweep_leading_edge_planform1=20,
#             wing_sweep_leading_edge_planform2=20,
#
#             front_spar_position=0.2,
#             rear_spar_position=0.6,
#             rib_number=12,
#
#             section_number=14,
#             points_number=14,
#             element_length=0.1
#         ), element_length=0.1,
#             thickness=thickness[0]
#         )
#
#         writer_new = Writer(instance_new)
#
#         # --- Paths ---
#         current_path = os.path.dirname(__file__)
#         output_dir = os.path.join(current_path, "output")
#         os.makedirs(output_dir, exist_ok=True)
#
#         mesh_path = os.path.join(output_dir, "mesh2.aster")
#         comm_path = os.path.join(output_dir, "wing_fem.comm")
#         export_path = os.path.join(output_dir, "case2.export")
#         results_path = os.path.join(output_dir, "results.txt")
#         log_file = os.path.join(output_dir, "aster_log.txt")
#
#         # --- Clean old files ---
#         for path in [mesh_path, comm_path, export_path, results_path, log_file]:
#             if os.path.exists(path):
#                 os.remove(path)
#
#         # --- Write new files ---
#         writer_new.write_mesh(mesh_path)
#         writer_new.write_comm(comm_path)
#         create_export_file(
#             target_path=export_path,
#             mesh_path=mesh_path,
#             comm_path=comm_path,
#             results_path=results_path,
#         )
#
#         # --- Run Code_Aster ---
#         run_code_aster(export_path, log_file=log_file)
#
#         # --- Read results ---
#         results_reader = ResultsReader(file_path=results_path)
#         max_deflection = float(results_reader.max_deflection)
#
#         print("Max deflection =", max_deflection)
#
#         # Constraint: max_deflection must be <= target_deflection
#         return target_deflection - max_deflection
#
#     # --- Define constraints for minimize ---
#     constraints = ({
#         'type': 'ineq',  # constraint >= 0
#         'fun': constraint
#     })
#
#     # --- Initial guess ---
#     x0 = [initial_thickness]
#
#     # --- Run the optimization ---
#     result = minimize(
#         objective,
#         x0,
#         method='SLSQP',  # good for constraints
#         bounds=[thickness_bounds],
#         constraints=constraints,
#         options={
#             'disp': True,
#             'ftol': 1e-9,
#             'maxiter': 100
#         }
#     )
#
#     return {
#         'optimized_thickness': result.x[0],
#         'max_deflection': target_deflection - constraint([result.x[0]]),
#         'success': result.success,
#         'message': result.message,
#         'result': result
#     }

def optimize_plate_thickness(target_deflection: float,
                             check_element,  # New argument for WingFEM instance
                             wing_airfoil_root,
                             wing_airfoil_middle,
                             wing_airfoil_tip,

                             wing_root_chord,
                             wing_middle_chord,
                             wing_tip_chord,

                             wing_thickness_factor_root,
                             wing_thickness_factor_middle,
                             wing_thickness_factor_tip,

                             wing_semi_span_planform1,
                             wing_semi_span,
                             wing_sweep_leading_edge_planform1,
                             wing_sweep_leading_edge_planform2,

                             front_spar_position,
                             rear_spar_position,
                             rib_number,

                             section_number,
                             points_number,
                             element_length,
                             avl_analysis,
                             skin_writer,
                             material_choice,
                             initial_thickness=0.1,
                             thickness_bounds=(0.02, 0.5),
                             ):
    """
    Optimize plate thickness to minimize thickness while keeping max deflection <= target_deflection.

    Args:
        target_deflection (float): Maximum allowed deflection (e.0.01 meters).
        wing_fem_instance: An instance of the WingFEM class.
        writer_instance: An instance of the Writer class.
        initial_thickness (float): Initial guess for the plate thickness.
        thickness_bounds (tuple): (min_thickness, max_thickness) allowed bounds.

    Returns:
        dict: Contains optimized thickness, max deflection, success status, and optimizer info.
    """

    def objective(thickness):
        """Objective: minimize thickness."""
        return thickness[0]  # minimize thickness itself

    def constraint(thickness):
        """Constraint: deflection must be <= target_deflection."""
        # --- Update WingFEM instance with new thickness ---
        # wing_fem_instance.thickness = thickness[0]

        print("Thickness =", thickness[0])

        instance_new = WingFEM(finalmesh=MeshGenerator(
            check_element=check_element,
            wing_airfoil_root=wing_airfoil_root,
            wing_airfoil_middle=wing_airfoil_middle,
            wing_airfoil_tip=wing_airfoil_tip,

            wing_root_chord=wing_root_chord,
            wing_middle_chord=wing_middle_chord,
            wing_tip_chord=wing_tip_chord,

            wing_thickness_factor_root=wing_thickness_factor_root,
            wing_thickness_factor_middle=wing_thickness_factor_middle,
            wing_thickness_factor_tip=wing_thickness_factor_tip,

            wing_semi_span_planform1=wing_semi_span_planform1,
            wing_semi_span=wing_semi_span,
            wing_sweep_leading_edge_planform1=wing_sweep_leading_edge_planform1,
            wing_sweep_leading_edge_planform2=wing_sweep_leading_edge_planform2,

            front_spar_position=front_spar_position,
            rear_spar_position=rear_spar_position,
            rib_number=rib_number,

            section_number=section_number,
            points_number=points_number,
            element_length=element_length
        ), element_length=element_length,
            thickness=thickness[0],
            avl=avl_analysis,
            skin_writer=skin_writer,
        )

        writer_new = Writer(instance=instance_new, avl=avl_analysis, material=material_choice)

        # --- Paths ---
        current_path = os.path.dirname(__file__)
        output_dir = os.path.join(current_path, "output")
        os.makedirs(output_dir, exist_ok=True)

        mesh_path = os.path.join(output_dir, "mesh2.aster")
        comm_path = os.path.join(output_dir, "wing_fem.comm")
        export_path = os.path.join(output_dir, "case2.export")
        results_path = os.path.join(output_dir, "results.txt")
        log_file = os.path.join(output_dir, "aster_log.txt")

        # --- Clean old files ---
        for path in [mesh_path, comm_path, export_path, results_path, log_file]:
            if os.path.exists(path):
                os.remove(path)

        # --- Write new files ---
        writer_new.write_mesh(mesh_path)
        writer_new.write_comm(comm_path)
        create_export_file(
            target_path=export_path,
            mesh_path=mesh_path,
            comm_path=comm_path,
            results_path=results_path,
        )

        # --- Run Code_Aster ---
        run_code_aster(export_path, log_file=log_file)

        # --- Read results ---
        results_reader = ResultsReader(file_path=results_path)
        max_deflection = float(results_reader.max_deflection)

        print("Max deflection =", max_deflection)

        # Constraint: max_deflection must be <= target_deflection
        return target_deflection - max_deflection

    # --- Define constraints for minimize ---
    constraints = ({
        'type': 'ineq',  # constraint >= 0
        'fun': constraint
    })

    # --- Initial guess ---
    x0 = [initial_thickness]

    # --- Run the optimization ---
    result = minimize(
        objective,
        x0,
        method='SLSQP',  # good for constraints
        bounds=[thickness_bounds],
        constraints=constraints,
        options={
            'disp': True,
            'ftol': 1e-9,
            'maxiter': 100
        }
    )

    return {
        'optimized_thickness': result.x[0],
        'max_deflection': target_deflection - constraint([result.x[0]]),
        'success': result.success,
        'message': result.message,
        'result': result
    }


if __name__ == "__main__":
    instance = WingFEM(finalmesh=MeshGenerator(
        check_element=0,
        wing_airfoil_root="whitcomb_interpolated.dat",
        wing_airfoil_middle="whitcomb_interpolated.dat",
        wing_airfoil_tip="whitcomb_interpolated.dat",

        wing_root_chord=6,
        wing_middle_chord=4,
        wing_tip_chord=1.5,

        wing_thickness_factor_root=1,
        wing_thickness_factor_middle=1,
        wing_thickness_factor_tip=1,

        wing_semi_span_planform1=5,
        wing_semi_span=16,
        wing_sweep_leading_edge_planform1=20,
        wing_sweep_leading_edge_planform2=20,

        front_spar_position=0.2,
        rear_spar_position=0.6,
        rib_number=12,

        section_number=14,
        points_number=14,
        element_length=0.1
    ), element_length=0.1
    )
    writer = Writer(instance)

    result = optimize_plate_thickness(target_deflection=0.07,
                                      wing_fem_instance=instance,
                                      writer_instance=writer)

    print(f"Optimized thickness: {result['optimized_thickness']} m")
    print(f"Max deflection achieved: {result['max_deflection']} m")
    print(f"Optimization success: {result['success']}")
