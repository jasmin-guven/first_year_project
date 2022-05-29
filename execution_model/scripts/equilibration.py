import BioSimSpace as bss
import sys
from BioSimSpace import _Exceptions
import os
from tqdm import tqdm

def run_process(system, md_protocol):
    """
    @param system: solvated system (bss object)
    @param md_protocol: bss protocolbss protocol
    @return the processed system.
    """
    process = bss.Process.Gromacs(system, md_protocol)
    process.setArg("-nt", 1)
    print(process.workDir())
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


minimisation_steps = 500
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = int(sys.argv[1])

ligand_stream = open("../ligands.dat", "r")
ligand_lines = ligand_stream.readlines()
ligand_name = ligand_lines[index].rstrip()

ligand_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_ligand_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_ligand_solvated.rst7"])
system_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_system_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_system_solvated.rst7"])
print("---------------------------")
print("Working on solvated ligand.")
print(f"Minimising in {minimisation_steps} steps.")
protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
for i in tqdm(range(minimisation_steps)):
    minimised = run_process(ligand_solvated, protocol)
print("Finished minimisation.")

print(f"NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = bss.Protocol.Equilibration(
                                    runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                    temperature_start=0*bss.Units.Temperature.kelvin,
                                    temperature_end=300*bss.Units.Temperature.kelvin,
                                    restraint="all"
                                )
for i in tqdm(range(runtime_short_nvt)):
    restrained_nvt = run_process(minimised, protocol)

print(f"NVT equilibration for {runtime_nvt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
for i in tqdm(range(runtime_nvt)):
    nvt = run_process(restrained_nvt, protocol)

print(f"NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )
for i in tqdm(range(runtime_npt)):
    restrained_npt = run_process(nvt, protocol)

print(f"NPT equilibration for {runtime_npt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )

for i in tqdm(range(runtime_npt)):
    equilibrated_ligand = run_process(restrained_npt, protocol)

print("---------------------------")
print("Working on solvated ligand+protein.")
print(f"Minimising in {minimisation_steps} steps..")
protocol = bss.Protocol.Minimisation(steps=minimisation_steps)

for i in tqdm(range(minimisation_steps)):
    minimised = run_process(system_solvated, protocol)

print(f"NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                temperature_start=0*bss.Units.Temperature.kelvin,
                                temperature_end=300*bss.Units.Temperature.kelvin,
                                restraint="all"
                                )

for i in tqdm(range(runtime_short_nvt)):
    restrained_nvt_system = run_process(minimised, protocol)

print(f"NVT equilibration for {runtime_nvt} ps while restraining all backbone atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="backbone"
                                )

for i in tqdm(range(runtime_nvt)):
    backbone_restrained_nvt_system = run_process(restrained_nvt_system, protocol)

print(f"NVT equilibration for {runtime_nvt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature_end=300*bss.Units.Temperature.kelvin,
                                )
for i in tqdm(range(runtime_nvt)):
    nvt_system = run_process(backbone_restrained_nvt_system, protocol)

print(f"NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms..")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )

for i in tqdm(range(runtime_npt)):
    restrained_npt_system = run_process(nvt_system, protocol)

print(f"NPT equilibration for {runtime_npt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
for i in tqdm(range(runtime_npt)):
    equilibrated_system = run_process(restrained_npt_system, protocol)

os.system("mkdir -p ../prep/ligands")
os.system("mkdir -p ../prep/protein")

print("Saving solvated/equilibrated systems.")
print("Ligand:")
print(equilibrated_ligand)
bss.IO.saveMolecules(f"../prep/ligands/{ligand_name}_ligand_prepped", equilibrated_ligand, ["PRM7", "RST7"])

print("\n Ligand + protein:")
print(equilibrated_system)
bss.IO.saveMolecules(f"../prep/protein/{ligand_name}_system_prepped", equilibrated_system, ["PRM7", "RST7"])
print("First 20 molecules in ligand + protein system:")
for molecules in equilibrated_system.getMolecules()[:20]:
    print(molecules)
print("Done.")
