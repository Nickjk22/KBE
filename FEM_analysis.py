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

# --- General imports ---
from collections import namedtuple
from typing import Dict, List, NamedTuple
from scipy.optimize import minimize
import os

# --- ParaPy core and geometry imports ---
from parapy.core import Attribute, Base, Input, Part, child
from glueing import GeneralFuse
from parapy.mesh import EdgeGroup, FaceGroup
from parapy.mesh.salome import Mesh
from parapy.geom import Compound, Face
from parapy.mesh.salome import MeshNode
from parapy.lib.code_aster import (_F, AFFE_CARA_ELEM, AFFE_CHAR_MECA,
                                   AFFE_MATERIAU, AFFE_MODELE, DEFI_MATERIAU,
                                   IMPR_RESU, LIRE_MAILLAGE, MECA_STATIQUE,
                                   Command, CommandWriter, MeshGroup,
                                   MeshWriter, ResultsReaderBase,
                                   create_export_file, run_code_aster)

# --- Project-specific modules ---
from meshing import MeshGenerator
from AVL_analysis import WingAVLAnalysis
from find_nodes import CodeAster_primitives
from wing import WingSurface

# --- Constants for constraint group labels ---
CONSTRAINED_EDGE1 = "constrained_edge1_group"
CONSTRAINED_EDGE2 = "constrained_edge2_group"
CONSTRAINED_EDGE3 = "constrained_edge3_group"
CONSTRAINED_EDGE4 = "constrained_edge4_group"


class WingFEM(Base):
    # Input parameters defining the FEM geometry and mesh
    thickness: float = Input(0.1)
    length: float = Input(1)
    width: float = Input(2)
    element_length: float = Input(0.1)

    @Input
    def finalmesh(self):
        # Mesh generator input
        return MeshGenerator(element_length=self.element_length)

    @Attribute
    def shape_to_mesh(self) -> GeneralFuse:
        # Geometry to mesh conversion
        return self.finalmesh.shape_to_mesh

    @Part
    def contrained_edge_groups(self) -> EdgeGroup:
        # Create edge groups used for constraints
        return EdgeGroup(quantify=4,
                         shape=[self.shape_to_mesh.edges[3],
                                self.shape_to_mesh.edges[6],
                                self.shape_to_mesh.edges[9],
                                self.shape_to_mesh.edges[10]][child.index],
                         label=[CONSTRAINED_EDGE1,
                                CONSTRAINED_EDGE2,
                                CONSTRAINED_EDGE3,
                                CONSTRAINED_EDGE4][child.index])

    @Attribute
    def FACE(self):
        # Generate labels for face groups
        return [f"face{i + 1}_group" for i in range(len(self.shape_to_mesh.faces))]

    @Part
    def face_group(self) -> FaceGroup:
        # Create face groups for shell definition
        return FaceGroup(quantify=len(self.shape_to_mesh.faces),
                         shape=[self.shape_to_mesh.faces[child.index]],
                         label=self.FACE[child.index])

    @Attribute
    def mesh(self) -> Mesh:
        # The generated mesh object
        return self.finalmesh.mesh

    @Input
    def avl(self):
        # AVL analysis input with basic setup
        return WingAVLAnalysis(aircraft=WingSurface(label="wing"), case_settings=[
            ("alpha_5deg", {'alpha': 5.0}),
        ])

    @Input
    def skin_writer(self):
        # Helper that provides load application nodes
        return CodeAster_primitives()

    @Attribute
    def liftforces(self) -> List[float]:
        # Extract lift forces from AVL
        return self.avl.lift_forces

    @Attribute
    def nodes(self) -> List[MeshNode]:
        # Load application nodes
        return self.skin_writer.load_primitives


class Writer:
    """Writes code_aster command and mesh files."""

    def __init__(self, instance: WingFEM = None, avl: WingAVLAnalysis = None, material: DEFI_MATERIAU = None) -> None:
        # If no WingFEM instance is passed, use default
        self._instance: WingFEM = instance or self._default_instance()
        # If no AVL is passed, use fallback setup
        self.avl: WingAVLAnalysis = avl or self._default_avl()
        # If no material is passed, use default aluminum
        self.material: DEFI_MATERIAU = material or self._default_material()

        # These will be filled by helper methods
        self.mesh_settings_command: LIRE_MAILLAGE = None
        self.mesh_groups: List[MeshGroup] = None
        self.model_command: AFFE_MODELE = None
        self.shell_properties_command: AFFE_CARA_ELEM = None
        self.material_zone_command: AFFE_MATERIAU = None
        self.constraint_commands: List[AFFE_CHAR_MECA] = None
        self.load_commands: List[AFFE_CHAR_MECA] = None
        self.solver_settings_command: MECA_STATIQUE = None
        self.result_settings_command: IMPR_RESU = None

        # Generate all Code_Aster setup commands
        self._generate_mesh_groups()
        self._generate_mesh_settings_command()
        self._generate_model_command()
        self._generate_shell_properties_command()
        self._generate_material_zone_command()
        self._generate_constraint_commands()
        self._generate_load_commands()
        self._generate_solver_settings_command()
        self._generate_results_settings_command()

        # Top-level commands for writing .comm
        self._commands: List[Command] = [self.result_settings_command]

    @staticmethod
    def _default_instance() -> WingFEM:
        # Creates default WingFEM for standalone execution
        return WingFEM()

    @staticmethod
    def _default_avl() -> WingAVLAnalysis:
        # Default AVL setup with simple alpha=5deg case
        wing = WingSurface(label="wing")
        return WingAVLAnalysis(aircraft=wing, case_settings=[("alpha_5deg", {'alpha': 5.0})])

    @staticmethod
    def _default_material() -> DEFI_MATERIAU:
        # Default material is aluminum
        ALUMINIUM = DEFI_MATERIAU(ELAS=_F(E=7e10, RHO=2700, NU=0.33))
        return ALUMINIUM

    @Attribute
    def liftforces(self):
        # Re-expose AVL lift forces
        return self.avl.lift_forces

    def write_comm(self, pathname: str) -> None:
        # Write the .comm file using Code_Aster CommandWriter
        writer = CommandWriter(self._commands)
        return writer.write(pathname)

    def write_mesh(self, pathname: str) -> None:
        # Write the .aster mesh file
        writer = MeshWriter(grid=self._instance.mesh.grid,
                            groups=self.mesh_groups)
        return writer.write(pathname)

    def _generate_mesh_groups(self) -> None:
        # Create MeshGroup objects for constraint edges and faces
        instance = self._instance

        # For each constrained edge, find corresponding mesh node group
        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE1,
                                                        shape=self._instance.shape_to_mesh.edges[3]).nodes
        ids = [elm.mesh_id for elm in elements]
        group1 = MeshGroup(label=CONSTRAINED_EDGE1, header="GROUP_NO", element_ids=ids, element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE2,
                                                        shape=self._instance.shape_to_mesh.edges[6]).nodes
        ids = [elm.mesh_id for elm in elements]
        group2 = MeshGroup(label=CONSTRAINED_EDGE2, header="GROUP_NO", element_ids=ids, element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE3,
                                                        shape=self._instance.shape_to_mesh.edges[9]).nodes
        ids = [elm.mesh_id for elm in elements]
        group3 = MeshGroup(label=CONSTRAINED_EDGE3, header="GROUP_NO", element_ids=ids, element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE4,
                                                        shape=self._instance.shape_to_mesh.edges[10]).nodes
        ids = [elm.mesh_id for elm in elements]
        group4 = MeshGroup(label=CONSTRAINED_EDGE4, header="GROUP_NO", element_ids=ids, element_type="node")

        # Store edge and face groups
        self.mesh_groups = [group1, group2, group3, group4]

        self.FACE = [f"face{i + 1}_group" for i in range(len(self._instance.shape_to_mesh.faces))]

        # Add face-based groups (surface elements)
        for i, face_label in enumerate(self.FACE):
            shape = self._instance.shape_to_mesh
            element_ids = []

            if isinstance(shape, Compound):
                faces = shape.faces
            elif isinstance(shape, Face):
                faces = [shape]
            else:
                faces = getattr(shape, 'faces', [])

            for face in faces:
                try:
                    subgrid = instance.mesh.get_subgrid_on_the_fly(label=face_label, shape=face)
                    element_ids.extend(f.mesh_id for f in subgrid.faces)
                except Exception as e:
                    print(f"Warning: Failed to process face for label {face_label}: {e}")

            group = MeshGroup(label=face_label, header="GROUP_MA", element_ids=element_ids, element_type="face")
            self.mesh_groups.append(group)

        # Add node groups for each load application point
        for i, mesh_node in enumerate(self._instance.nodes):
            grp_label = f"load_node{i + 1}_group"
            node_group = MeshGroup(
                label=grp_label,
                header="GROUP_NO",
                element_ids=[mesh_node.mesh_id],
                element_type="node"
            )
            self.mesh_groups.append(node_group)

    def _generate_mesh_settings_command(self) -> None:
        # Code_Aster command for mesh settings
        self.mesh_settings_command = LIRE_MAILLAGE(UNITE=20, FORMAT="ASTER")

    def _generate_model_command(self) -> None:
        # Assign DKT shell elements to faces
        shell_model = _F(GROUP_MA=tuple(self.FACE),
                         MODELISATION=("DKT",),
                         PHENOMENE="MECANIQUE")
        self.model_command = AFFE_MODELE(AFFE=(shell_model,),
                                         MAILLAGE=self.mesh_settings_command)

    def _generate_shell_properties_command(self) -> None:
        # Assign thickness and local direction vector to shell elements
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
        # Assign material properties to shell elements
        self.material_zone_command = AFFE_MATERIAU(AFFE=(_F(GROUP_MA=tuple(self.FACE),
                                                            MATER=(self.material,)),),
                                                   MAILLAGE=self.mesh_settings_command,
                                                   MODELE=self.model_command)

    def _generate_constraint_commands(self) -> None:
        # Clamp the four edge groups with ENCASTRE boundary condition
        self.constraint_commands = [AFFE_CHAR_MECA(
            DDL_IMPO=_F(
                GROUP_NO=(CONSTRAINED_EDGE1, CONSTRAINED_EDGE2, CONSTRAINED_EDGE3, CONSTRAINED_EDGE4),
                LIAISON="ENCASTRE"),
            MODELE=self.model_command)]

    def _generate_load_commands(self) -> None:
        # Apply nodal forces (lift) to each load node group
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
        # Create MECA_STATIQUE solver command combining loads and constraints
        loads_constraints = [_F(CHARGE=obj) for obj in self.load_commands + self.constraint_commands]
        self.solver_settings_command = MECA_STATIQUE(SOLVEUR=_F(RESI_RELA=1e-06),
                                                     CARA_ELEM=self.shell_properties_command,
                                                     CHAM_MATER=self.material_zone_command,
                                                     EXCIT=loads_constraints,
                                                     MODELE=self.model_command)

    def _generate_results_settings_command(self):
        # Specify output of displacement results
        self.result_settings_command = IMPR_RESU(FORMAT="RESULTAT",
                                                 RESU=(_F(NOM_CHAM=("DEPL",),
                                                          RESULTAT=self.solver_settings_command),),
                                                 UNITE=8)


class ResultsReader(ResultsReaderBase):
    @Attribute
    def nodes_displacements(self) -> Dict[str, NamedTuple]:
        # Parse displacement results into a dictionary of namedtuples
        results_ = self.results["DEPL"]
        tmplt = namedtuple("nodes_displacements", ["TX", "TY", "TZ", "RX", "RY", "RZ"])
        dct = {}
        for line in results_[1:]:
            line_ = line.split()
            dct[line_[0]] = tmplt(*line_[1:])
        return dct

    @Attribute
    def max_deflection(self) -> float:
        # Compute the maximum deflection in the Z direction
        deflections = [float(value.TZ) for value in self.nodes_displacements.values()]
        return max(deflections)


def optimize_plate_thickness(target_deflection: float,
                             check_element,
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
    Optimize plate thickness to minimize weight while keeping deflection under the target.

    Args:
        target_deflection (float): Upper limit on allowed max deflection (e.g., 0.01 meters).
        check_element (int): Identifier for meshing config.
        wing_airfoil_*: Airfoil definitions for root, middle, and tip.
        wing_*: Geometry and configuration parameters of the wing.
        avl_analysis: AVL force results used as loads.
        skin_writer: Provides nodes to apply loads.
        material_choice: Material to be used in the model.
        initial_thickness (float): Starting point for optimization.
        thickness_bounds (tuple): Min and max bounds on thickness.

    Returns:
        dict: Includes optimized thickness, deflection, success flag, and full result object.
    """

    def objective(thickness):
        # Objective: minimize the thickness itself
        return thickness[0]

    def constraint(thickness):
        # Constraint: max deflection should not exceed target
        print("Thickness =", thickness[0])

        # Create a new WingFEM instance with updated thickness
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

        # --- Define file paths ---
        current_path = os.path.dirname(__file__)
        output_dir = os.path.join(current_path, "output")
        os.makedirs(output_dir, exist_ok=True)

        mesh_path = os.path.join(output_dir, "mesh2.aster")
        comm_path = os.path.join(output_dir, "wing_fem.comm")
        export_path = os.path.join(output_dir, "case2.export")
        results_path = os.path.join(output_dir, "results.txt")
        log_file = os.path.join(output_dir, "aster_log.txt")

        # --- Clean up old files if they exist ---
        for path in [mesh_path, comm_path, export_path, results_path, log_file]:
            if os.path.exists(path):
                os.remove(path)

        # --- Generate new files for Code_Aster ---
        writer_new.write_mesh(mesh_path)
        writer_new.write_comm(comm_path)
        create_export_file(
            target_path=export_path,
            mesh_path=mesh_path,
            comm_path=comm_path,
            results_path=results_path,
        )

        # --- Run simulation ---
        run_code_aster(export_path, log_file=log_file)

        # --- Read results and extract max deflection ---
        results_reader = ResultsReader(file_path=results_path)
        max_deflection = float(results_reader.max_deflection)

        print("Max deflection =", max_deflection)

        # Constraint: target - max_deflection >= 0
        return target_deflection - max_deflection

    # Define constraint for optimizer
    constraints = ({
        'type': 'ineq',  # Must be >= 0
        'fun': constraint
    })

    # Starting value for thickness
    x0 = [initial_thickness]

    # Run the optimization using Sequential Least Squares
    result = minimize(
        objective,
        x0,
        method='SLSQP',
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
    # Example usage of the WingFEM and optimization
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
