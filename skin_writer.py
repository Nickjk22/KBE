from skin import CodeAster_primitives
from parapy.lib.code_aster import (_F, AFFE_CARA_ELEM, AFFE_CHAR_MECA,
                                   AFFE_MATERIAU, AFFE_MODELE, DEFI_MATERIAU,
                                   IMPR_RESU, LIRE_MAILLAGE, MECA_STATIQUE,
                                   Command, CommandWriter, MeshGroup,
                                   MeshWriter, ResultsReaderBase,
                                   create_export_file, run_code_aster)

from os.path import join
import os

class SkinWriter:
    def __init__(self, instance: CodeAster_primitives, working_dir="output/skin_run"):
        self.instance = instance
        self.working_dir = working_dir
        os.makedirs(working_dir, exist_ok=True)

        self.mesh = self.instance.finalmesh.mesh_generator.mesh.grid
        self.subgrids = self.instance.subgrids
        self.primitives = self.instance.primitives
        self.load_primitives = self.instance.load_primitives

        self.commands = []

    def write(self):
        mesh_path = join(self.working_dir, "mesh.med")
        comm_path = join(self.working_dir, "study.comm")
        export_path = join(self.working_dir, "study.export")
        log_path = join(self.working_dir, "study.log")

        # Write mesh
        load_nodes = [node.id for node in self.load_primitives]
        load_group = MeshGroup("load_nodes", load_nodes, "NOEUD")
        MeshWriter(self.mesh, groups=self.subgrids + [load_group]).write(mesh_path)

        # Build commands
        model = AFFE_MODELE(
            mail=self.mesh,
            groupma=self.subgrids,
            model="3D",
            carac="MECANIQUE"
        )

        material = DEFI_MATERIAU(elas=Elas(E=70000e6, NU=0.33))

        fieldmat = AFFE_MATERIAU(
            mail=self.mesh,
            model=model,
            mater=material
        )

        meshfunc = AFFE_CARA_ELEM(
            mail=self.mesh,
            model=model,
            penta=PENTA(epa=0.005)
        )

        loads = [
            LOAD(
                groupno=load_prim,
                fx=0,
                fy=0,
                fz=-1000
            )
            for load_prim in self.load_primitives
        ]

        force = AFFE_CHAR_MECA(
            model=model,
            force=loads
        )

        statics = MECA_STATIQUE(
            model=model,
            cara_elem=meshfunc,
            mater=fieldmat,
            char=[force]
        )

        result = IMPR_RESU(
            FORMAT="MED",
            RESU=[
                Result(_resultat=statics, nom_champ="DEPL", NOM_CHAM="DEPL", RESULTAT=statics)
            ]
        )

        self.commands = [model, material, fieldmat, meshfunc, force, statics, result]

        # Write .comm file
        CommandWriter(self.commands).write(comm_path)

        # Write .export file
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
