# solvation.py
import BioSimSpace as bss
from BioSimSpace import _Exceptions
import sys
import os
import csv

minimisation_steps = 250
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = int(sys.argv[1])

ligand_stream = open("../ligands.dat", "r")
ligand_lines = ligand_stream.readlines()
ligand_name = ligand_lines[index].rstrip()

ligand = bss.IO.readMolecules(f"inputs/ligands/{ligand_name}.mol2")[0]

protocol_stream = open("../protocol.dat", "r")
protocol_lines = protocol_stream.readlines()
ligand_force_field = protocol_lines[0].rstrip()

if "GAFF2" in ligand_force_field:
    print(f"Parameterising {ligand_nameg} using GAFF2 force field.")
    ligand_params = BSS.Parameters.gaff2(ligand).getMolecule()
else:
    raise NameError(f"Force field not supported: {ligand_force_field}. Please use either of [GAFF1, GAFF2, "
                    f"OpenForceField], or edit this script to use other force fields available in "
                    f"BSS.Parameters.forceFields().")
ligand_params_copy = ligand_params.copy()

solvent_force_field = protocol_lines[2].rstrip().replace(" ", "").split("=")[-1]
box_size = protocol_lines[3].rstrip().replace(" ", "").split("=")[-1]
box_axis_length = box_size.split("*")[0]
box_axis_unit = BSS.Units.Length.angstrom
box_type = protocol_lines[4].rstrip().replace(" ", "").split("=")[-1]
print(f"box type is set to {box_type}")
"""
figure out the ligand dimensions, add the specified water box together with the largest ligand axis.
This is a workaround to account for adding ions in some cases. Based on:
https://github.com/michellab/BioSimSpace/issues/111
"""
box_min, box_max = ligand_params.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]

protein = BSS.IO.readMolecules(["../inputs/protein/protein_complex.rst7", "../inputs/protein/protein_complex.prm7"])[0]
system = ligand_params + protein
box_min_system, box_max_system = system.getAxisAlignedBoundingBox()
box_size_system = [y - x for x, y in zip(box_min_s, box_max_system)]
box_sizes_system = [x + int(box_axis_length) * box_axis_unit for x in box_size_system]

print("Solvating ligand.")
box, angles = BSS.Box.cubic(max(box_sizes))
ligand_params_solvated = BSS.Solvent.solvate(solvent_force_field, molecule=ligand_params, box=box, angles=angles)

print("Solvating ligand + protein system.")
box, angles = BSS.Box.cubic(max(box_sizes_s))
system_solvated = BSS.Solvent.solvate(solvent_force_field, molecule=system, box=box, angles=angles)

BSS.IO.saveMolecules(f"temp/{ligand_name}_ligand_solvated", ligand_params_solvated, ["PRM7", "RST7"])
BSS.IO.saveMolecules(f"temp/{ligand_name}_system_solvated", system_solvated, ["PRM7", "RST7"])
