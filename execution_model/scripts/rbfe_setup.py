# Relative Binding free-energy (RBFE) calculation setup
import BioSimSpace as bss
from BioSimSpace import _Exceptions
import os
import sys
import csv
import warnings
warnings.filterwarnings("ignore")

print(f"Program: {sys.argv[0]}, ligand 1: {sys.argv[1]}, ligand 2: {sys.argv[2]}")

print(f"Loading ligands {sys.argv[1]} and {sys.argv[2]}.")
ligand_path = "../prep/ligands/"
#index = int(sys.argv[1]) + 1 
#print(f"index: {index}")
ligand_1_sys = bss.IO.readMolecules([f"{ligand_path}lig_h_{sys.argv[1]}_ligand_prepped.rst7",
                                     f"{ligand_path}lig_h_{sys.argv[1]}_ligand_prepped.prm7"])
ligand_2_sys = bss.IO.readMolecules([f"{ligand_path}lig_h_{sys.argv[2]}_ligand_prepped.rst7",
                                     f"{ligand_path}lig_h_{sys.argv[2]}_ligand_prepped.prm7"])
print(f"ligand_1_sys: {ligand_path}lig_h_{sys.argv[1]}_ligand_prepped.prm7")
print(f"ligand_2_sys: {ligand_path}lig_h_{sys.argv[2]}_ligand_prepped.prm7")
ligand_1 = ligand_1_sys.getMolecule(0)
ligand_2 = ligand_2_sys.getMolecule(0)

print("Mapping and aligning...")
print(ligand_1, ligand_2)
mapping = bss.Align.matchAtoms(ligand_1, ligand_2,  complete_rings_only=True)
inverse_mapping = {v: k for k, v in mapping.items()}
ligand_2_aligned = bss.Align.rmsdAlign(ligand_2, ligand_1, inverse_mapping)

print("Merging...")
try:
    merged_ligands = bss.Align.merge(ligand_1, ligand_2_aligned, mapping)
except _Exceptions.IncompatibleError:
    raise NameError(f"The merge has opened/closed a ring on ligand {sys.argv[2]}")
#### Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
# Following is work around because setBox() doesn't validate correctly boxes with lengths and angles

ligand_1_sys.removeMolecules(ligand_1)
ligand_1_sys.addMolecules(merged_ligands)
system_free = ligand_1_sys
print(f"system free: {system_free}")

# now repeat above steps, but for the protein + ligand systems.
# Load equilibrated bound inputs for both ligands. Complain if input not found
print(f"Loading bound ligands {sys.argv[1]} and {sys.argv[2]}.")
protein_path = "../prep/protein/"
system_1 = bss.IO.readMolecules(
    [f"{protein_path}lig_h_{sys.argv[1]}_system_prepped.rst7",
     f"{protein_path}lig_h_{sys.argv[1]}_system_prepped.prm7"])
system_2 = bss.IO.readMolecules(
    [f"{protein_path}lig_h_{sys.argv[2]}_system_prepped.rst7",
     f"{protein_path}lig_h_{sys.argv[2]}_system_prepped.prm7"])


# Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
# the order of molecules is switched, so we can't use index alone.
# bugfix in BSS makes the below redundant but keeping this in to be 100% sure we're getting the correct structures.

system_ligand_1 = None
protein = None

n_residues = [molecule.nResidues() for molecule in system_1]
n_atoms = [molecule.nAtoms() for molecule in system_1]

for i, (n_resi, n_at) in enumerate(zip(n_residues[:20], n_atoms[:20])):
    if n_resi == 1 and n_at > 5:
        # grab the system's ligand and the protein. ignore the waters.
        system_ligand_1 = system_1.getMolecule(i)
    elif n_resi > 1:
        protein = system_1.getMolecule(i)
    else:
        pass

# loop over molecules in system to extract the ligand
system_ligand_2 = None

print(f"system_ligand_1: {protein_path}lig_h_{sys.argv[1]}_system_prepped.rst7")
print(f"system_ligand_2: {protein_path}lig_h_{sys.argv[2]}_system_prepped.rst7")

n_residues = [mol.nResidues() for mol in system_2]
n_atoms = [mol.nAtoms() for mol in system_2]
print(n_residues)
for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):

    # grab the system's ligand. ignore the waters.
    if n_resi == 1 and n_at > 5:
        system_ligand_2 = system_2.getMolecule(i)
    else:
        pass

if system_ligand_2:
    print("================")
    print("system ligand 2 OK")
    print("================")
    
if system_ligand_1 and system_ligand_2 and protein:
    print("Using molecules ligand_1, ligand_2, protein:")
    print(system_ligand_1, system_ligand_2, protein)
else:
    raise _Exceptions.AlignmentError(
        "Could not extract ligands or protein from input systems. "
        "Check that your ligands/proteins are properly prepared by BSSligprep.sh!")

# Align ligand2 on ligand1
print("Mapping...")
mapping = bss.Align.matchAtoms(system_ligand_1, system_ligand_2, complete_rings_only=True)
inverse_mapping = {v: k for k, v in mapping.items()}

print("Aligning...")
system_ligand_2_aligned = bss.Align.rmsdAlign(system_ligand_2, system_ligand_1, inverse_mapping)

# Generate merged molecule.
print("Merging...")
system_merged_ligands = bss.Align.merge(system_ligand_1, system_ligand_2_aligned, mapping)

system_1.removeMolecules(system_ligand_1)
system_1.addMolecules(system_merged_ligands)
system_bound = system_1

stream = open("../protocol.dat", "r")
lines = stream.readlines()

engine = lines[7].rstrip().replace(" ", "").split("=")[-1].upper()
if engine not in ["SOMD", "GROMACS"]:
    raise NameError("Input MD engine not recognised. Please use any of ['SOMD', 'GROMACS']" \
                    + "on the eighth line of protocol.dat in the shape of (e.g.):\nengine = SOMD")

runtime = lines[6].rstrip().replace(" ", "").split("=")[-1].split("*")[0]
try:
    runtime = int(runtime)
except ValueError:
    raise NameError("Input runtime value not supported. Please use an integer" \
                    + " on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

# make sure user has set ns or ps.
runtime_unit = lines[6].rstrip().replace(" ", "").split("=")[-1].split("*")[1]
if runtime_unit not in ["ns", "ps"]:
    raise NameError("Input runtime unit not supported. Please use 'ns' or 'ps'" \
                    + " on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")

if runtime_unit == "ns":
    runtime_unit = bss.Units.Time.nanosecond
elif runtime_unit == "ps":
    runtime_unit = bss.Units.Time.picosecond

# get the number of lambda windows for this pert.
n_windows = None
ligand_1_name = f"lig_h_{sys.argv[1]}"
ligand_2_name = f"lig_h_{sys.argv[2]}"
with open("../network.dat", "r") as lambdas_file:
    reader = csv.reader(lambdas_file, delimiter=" ")
    for row in reader:
        if row[0] == ligand_1_name and row[1] == ligand_2_name:
            n_windows = int(row[2])
if not n_windows:
    raise NameError(f"The perturbation {ligand_1_name}~{ligand_2_name} was not found in network.dat.")

# define the free energy protocol with all this information. User could customise settings further here, see docs.
rbfe_protocol = bss.Protocol.FreeEnergy(num_lam=n_windows, runtime=runtime * runtime_unit)

# Set up the directory environment.
# testing is already done by BSS.
workdir = f"../outputs/{engine}/{ligand_1_name}~{ligand_2_name}/"
print(f"Setting up {engine} directory environment in {workdir}.")

# set up a bound folder with standard settings.
# Use solvation to prep the bound leg
print("Bound...")
workdir = f"../outputs/{engine}/{ligand_1_name}~{ligand_2_name}"
bss.FreeEnergy.Relative(
    system_bound,
    rbfe_protocol,
    engine=f"{engine}",
    work_dir=workdir + "/bound"
)

# set up a free folder.
print("Free...")
bss.FreeEnergy.Relative(
    system_free,
    rbfe_protocol,
    engine=f"{engine}",
    work_dir=workdir + "/free"
)
