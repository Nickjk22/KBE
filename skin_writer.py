from find_nodes import CodeAster_primitives
from parapy.lib.code_aster import (
    _F, AFFE_CARA_ELEM, AFFE_CHAR_MECA, AFFE_MATERIAU, AFFE_MODELE,
    DEFI_MATERIAU, IMPR_RESU, LIRE_MAILLAGE, MECA_STATIQUE,
    MeshWriter, CommandWriter, create_export_file, run_code_aster
)
from parapy.lib.code_aster import MeshGroup
from parapy.core import *

from os.path import join
import os


class SkinWriter:

    def __init__(self, instance: CodeAster_primitives, working_dir="output/skin_run"):
        self.instance = instance
        self.working_dir = working_dir
        os.makedirs(working_dir, exist_ok=True)

    def write(self):
        mesh = self.instance.finalmesh.mesh_generator.mesh
        mesh_path = join(self.working_dir, "mesh.med")
        comm_path = join(self.working_dir, "study.comm")
        export_path = join(self.working_dir, "study.export")
        log_path = join(self.working_dir, "study.log")

        # Generate subgrids from the CAD-to-mesh history
        # subgrids = self.instance.subgrids

        print(dir(self.instance.finalmesh.shape_to_mesh))

        shape = self.instance.finalmesh.shape_to_mesh.get_slot_value("TopoDS_Shape")
        print("Shape:", shape)

        if shape is None:
            raise ValueError("Could not retrieve valid shape_in from shape_to_mesh")

        history = self.instance.finalmesh.shape_to_mesh.history(shape)
        shapes = list(history.keys())

        # Zet de subgrids om naar MeshGroups (nodig voor MeshWriter)
        groups = [
            MeshGroup(name=sg.label, elements=sg.elements, nodes=sg.nodes)
            for sg in subgrids
        ]

        # Schrijf mesh met groepen naar bestand
        MeshWriter(mesh.grid, groups=groups).write(mesh_path)

        # Modeldefinitie
        model = AFFE_MODELE(
            MAILLAGE=mesh,
            AFFE=_F(GROUP_MA=[group.name for group in subgrids], MODELISATION="3D")
        )

        # Materiaaldefinitie
        material = DEFI_MATERIAU(ELAS=_F(E=70e9, NU=0.33))

        fieldmat = AFFE_MATERIAU(
            mail=mesh,
            model=model,
            mater=_F(GROUP_MA=subgrids[0].name, MATER=material)
        )

        # Dikte toekennen aan elke subgrid
        cara_elem = AFFE_CARA_ELEM(
            mail=mesh,
            model=model,
            penta=[_F(GROUP_MA=group.name, EPA=self.instance.skin.thickness)
                   for group in subgrids]
        )

        # Externe krachten
        loads = []
        for node, lift in zip(self.instance.load_primitives, self.instance.lift_forces):
            loads.append(
                _F(GROUP_NO=node.name, FX=0, FY=0, FZ=-lift)
            )

        force = AFFE_CHAR_MECA(
            model=model,
            force=loads
        )

        # Mechanische statische berekening
        statics = MECA_STATIQUE(
            model=model,
            cara_elem=cara_elem,
            mater=fieldmat,
            char=[force]
        )

        # Resultaat wegschrijven
        result = IMPR_RESU(
            FORMAT="MED",
            RESU=[
                _F(RESULTAT=statics, NOM_CHAM="DEPL", NOM_CMP=("DX", "DY", "DZ"))
            ]
        )

        # Commands schrijven
        commands = [model, material, fieldmat, cara_elem, force, statics, result]
        CommandWriter(commands).write(comm_path)

        # Exportfile schrijven
        create_export_file(
            comm_file=comm_path,
            mesh_file=mesh_path,
            result_dir=self.working_dir,
            export_file=export_path
        )

        return export_path, log_path

    def run(self):
        export_file, log_file = self.write()
        run_code_aster(export_file, log_file=log_file)
        print(f"Simulation complete. Log: {log_file}")


if __name__ == "__main__":
    skin_model = CodeAster_primitives()
    writer = SkinWriter(skin_model)
    writer.run()
