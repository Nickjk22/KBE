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
from aircraft_fem.examples.aircraft.geom.glueing import GeneralFuse
from parapy.geom import RectangularFace
from parapy.mesh import EdgeGroup, FaceGroup
from parapy.mesh.salome import Mesh, Tri
from parapy.geom import Compound, Face
from parapy.mesh.salome.controls import FixedLength
from parapy.lib.code_aster import (_F, AFFE_CARA_ELEM, AFFE_CHAR_MECA,
                                   AFFE_MATERIAU, AFFE_MODELE, DEFI_MATERIAU,
                                   IMPR_RESU, LIRE_MAILLAGE, MECA_STATIQUE,
                                   Command, CommandWriter, MeshGroup,
                                   MeshWriter, ResultsReaderBase,
                                   create_export_file, run_code_aster)
from meshing_riks import FinalMesh, MeshGenerator
from AVL_main import WingAVLAnalysis
from skin import CodeAster_primitives

# FACE = "face_group"
CONSTRAINED_EDGE1 = "constrained_edge1_group"
CONSTRAINED_EDGE2 = "constrained_edge2_group"
CONSTRAINED_EDGE3 = "constrained_edge3_group"
CONSTRAINED_EDGE4 = "constrained_edge4_group"
LOADED_EDGE = "loaded_edge_group"

STEEL = DEFI_MATERIAU(ELAS=_F(E=300000000000.0, RHO=7850, NU=0.1666))


class WingFEM(Base):
    element_length: float = Input(0.1)
    thickness: float = Input(0.1)
    length: float = Input(1)
    width: float = Input(2)

    # @Part
    # def shape_to_mesh(self) -> RectangularFace:
    #     return RectangularFace(length=self.length,
    #                            width=self.width)

    # @Part
    # def hyp_1d(self) -> FixedLength:
    #     return FixedLength(length=self.element_length,
    #                        shape=self.shape_to_mesh)
    #
    # @Part
    # def hyp_2d(self) -> Tri:
    #     return Tri(shape=self.shape_to_mesh)

    # @Part
    # def contrained_edge_groups(self) -> EdgeGroup:
    #     return EdgeGroup(quantify=2,
    #                      shape=self.shape_to_mesh.edges[0].neighbours[child.index],
    #                      label=[CONSTRAINED_EDGE1, CONSTRAINED_EDGE2][child.index])

    # @Part
    # def loaded_edge_group(self) -> EdgeGroup:
    #     return EdgeGroup(shape=self.shape_to_mesh.edges[0],
    #                      label=LOADED_EDGE)

    # @Part
    # def face_group(self) -> FaceGroup:
    #     return FaceGroup(shape=self.shape_to_mesh,
    #                      label=FACE)

    # @Part
    # def mesh(self) -> Mesh:
    #     return Mesh(shape_to_mesh=self.shape_to_mesh,
    #                 groups=[self.face_group,
    #                         self.loaded_edge_group] + self.contrained_edge_groups,
    #                 controls=[self.hyp_1d, self.hyp_2d])

    @Part
    def finalmesh(self):
        return FinalMesh()

    @Attribute
    def shape_to_mesh(self) -> GeneralFuse:
        return self.finalmesh.shape_to_mesh

    @Part
    def contrained_edge_groups(self) -> EdgeGroup:
        return EdgeGroup(quantify=4,
                         shape=[self.shape_to_mesh.edges[3], self.shape_to_mesh.edges[6], self.shape_to_mesh.edges[9],
                                self.shape_to_mesh.edges[10]][child.index],
                         label=[CONSTRAINED_EDGE1, CONSTRAINED_EDGE2, CONSTRAINED_EDGE3, CONSTRAINED_EDGE4][
                             child.index])

    @Part
    def loaded_edge_group(self) -> EdgeGroup:
        return EdgeGroup(shape=self.shape_to_mesh.edges[157],
                         label=LOADED_EDGE)

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

    def __init__(self, instance: WingFEM) -> None:
        self._instance: WingFEM = instance
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

    @Part
    def avl(self):
        return WingAVLAnalysis()

    @Part
    def skin_writer(self):
        return CodeAster_primitives()

    @Attribute
    def lift_forces(self):
        return self.avl.lift_forces

    @Attribute
    def nodes(self):
        return self.skin_writer.load_primitives

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
        ids = [elm.mesh_id for elm in elements if elm.y < instance.length / 2]
        group1 = MeshGroup(label=CONSTRAINED_EDGE1,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE2,
                                                        shape=self._instance.shape_to_mesh.edges[6]).nodes
        ids = [elm.mesh_id for elm in elements if elm.y < instance.length / 2]
        group2 = MeshGroup(label=CONSTRAINED_EDGE2,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE3,
                                                        shape=self._instance.shape_to_mesh.edges[9]).nodes
        ids = [elm.mesh_id for elm in elements if elm.y < instance.length / 2]
        group3 = MeshGroup(label=CONSTRAINED_EDGE3,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=CONSTRAINED_EDGE4,
                                                        shape=self._instance.shape_to_mesh.edges[10]).nodes
        ids = [elm.mesh_id for elm in elements if elm.y < instance.length / 2]
        group4 = MeshGroup(label=CONSTRAINED_EDGE4,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        elements = instance.mesh.get_subgrid_on_the_fly(label=LOADED_EDGE,
                                                        shape=self._instance.shape_to_mesh.edges[157]).nodes
        ids = [elm.mesh_id for elm in elements]
        group5 = MeshGroup(label=LOADED_EDGE,
                           header="GROUP_NO",
                           element_ids=ids,
                           element_type="node")

        # group6 = MeshGroup(label=FACE,
        #                    header="GROUP_MA",
        #                    element_ids=[f.mesh_id for f in instance.mesh.get_subgrid_on_the_fly(label=FACE,
        #                                                                                         shape=self._instance.shape_to_mesh).faces],
        #                    element_type="face")

        self.mesh_groups = [group1, group2, group3, group4, group5]

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
                                                            MATER=(STEEL,)),),
                                                   MAILLAGE=self.mesh_settings_command,
                                                   MODELE=self.model_command)

    def _generate_constraint_commands(self) -> None:
        self.constraint_commands = [AFFE_CHAR_MECA(
            DDL_IMPO=_F(
                GROUP_NO=(CONSTRAINED_EDGE1, CONSTRAINED_EDGE2, CONSTRAINED_EDGE3, CONSTRAINED_EDGE4),
                LIAISON="ENCASTRE"),
            MODELE=self.model_command)]

    def _generate_load_commands(self) -> None:
        self.load_commands = [AFFE_CHAR_MECA(FORCE_NODALE=_F(FZ=10,
                                                             GROUP_NO=(LOADED_EDGE,)),
                                             MODELE=self.model_command), AFFE_CHAR_MECA(FORCE_NODALE=_F(FZ=10,
                                                                                                        GROUP_NO=(
                                                                                                        LOADED_EDGE,)),
                                                                                        MODELE=self.model_command)]

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


if __name__ == "__main__":
    import os

    instance = WingFEM()
    # display(instance)

    current_path = os.path.dirname(__file__)
    output_dir = os.path.join(current_path, "output")

    writer = Writer(instance)
    mesh_path = os.path.join(output_dir, "mesh2.aster")
    writer.write_mesh(mesh_path)

    comm_path = os.path.join(output_dir, "wing_fem.comm")
    writer.write_comm(comm_path)

    export_path = os.path.join(output_dir, "case2.export")
    results_path = os.path.join(output_dir, "results.txt")
    create_export_file(
        target_path=export_path,
        mesh_path=mesh_path,
        comm_path=comm_path,
        results_path=results_path,
    )

    log_file = os.path.join(output_dir, "aster_log.txt")
    run_code_aster(export_path, log_file=log_file)

    results_reader = ResultsReader(file_path=results_path)
    print(f"maximum deflection from FEM: {float(results_reader.max_deflection)}")

#
# def optimize_plate_thickness(target_deflection: float, thickness_bounds=(0.01, 0.5)):
#     """
#     Optimize plate thickness to minimize thickness while keeping max deflection <= target_deflection.
#
#     Args:
#         target_deflection (float): Maximum allowed deflection (e.g., 0.01 meters).
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
#         # --- Set up new Plate ---
#         instance = Plate(thickness=thickness[0])
#
#         # --- Set up writer ---
#         writer = Writer(instance)
#
#         # --- Paths ---
#         current_path = os.path.dirname(__file__)
#         output_dir = os.path.join(current_path, "output")
#         os.makedirs(output_dir, exist_ok=True)
#
#         mesh_path = os.path.join(output_dir, "mesh2.aster")
#         comm_path = os.path.join(output_dir, "rectangular_plate.comm")
#         export_path = os.path.join(output_dir, "case2.export")
#         results_path = os.path.join(output_dir, "results.txt")
#         log_file = os.path.join(output_dir, "aster_log.txt")
#
#         # --- Write new files ---
#         writer.write_mesh(mesh_path)
#         writer.write_comm(comm_path)
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
#     x0 = [0.1]  # start at 0.1 meters thickness
#
#     # --- Run the optimization ---
#     result = minimize(
#         objective,
#         x0,
#         method='SLSQP',  # good for constraints
#         bounds=[thickness_bounds],
#         constraints=constraints,
#         options={'disp': True}
#     )
#
#     return {
#         'optimized_thickness': result.x[0],
#         'max_deflection': target_deflection - constraint([result.x[0]]),
#         'success': result.success,
#         'message': result.message,
#         'result': result
#     }
#
#
# if __name__ == "plate":
#     # Example: target deflection 0.01 meters (10 mm)
#     result = optimize_plate_thickness(target_deflection=0.0001)
#
#     print(f"Optimized thickness: {result['optimized_thickness']} m")
#     print(f"Max deflection achieved: {result['max_deflection']} m")
#     print(f"Optimization success: {result['success']}")
