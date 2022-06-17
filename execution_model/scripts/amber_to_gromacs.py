import BioSimSpace as bss
import sys
from BioSimSpace import _Exceptions
import os

def amb2gmx(system, savename):
    try:
        bss.IO.saveMolecules(savename, system, ["Gro87", "GroTop"])
    except:
        print("Could not convert.")


print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = int(sys.argv[1])

ligand_stream = open("../ligands.dat", "r")
ligand_lines = ligand_stream.readlines()
ligand_name = ligand_lines[index].rstrip()

ligand_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_ligand_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_ligand_solvated.rst7"])
system_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_system_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_system_solvated.rst7"])

amb2gmx(ligand_solvated, f"../inputs/ligands/{ligand_name}_ligand_solvated")
amb2gmx(system_solvated, f"../inputs/ligands/{ligand_name}_system_solvated")