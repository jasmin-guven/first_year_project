# equilibration.py

import BioSimSpace as BSS
import sys
import os
import csv

minimisation_steps = 250
runtime_short_nvt = 5 # ps
runtime_nvt = 50 # ps
runtime_npt = 200 # ps

print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = sys.argv[1]

file_stream = open("../execution_model/ligands.dat", "r")
lines = file_stream.readlines()
ligand_name = lines[index].rstrip()

print(ligand_name)