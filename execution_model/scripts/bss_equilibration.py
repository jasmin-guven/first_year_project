import BioSimSpace as bss
import sys
from BioSimSpace import _Exceptions
import os


def run_process(system, md_protocol):
    process = bss.Process.Gromacs(system, md_protocol, "/opt/bin/bin/gmx")
    process.setArg("-nt", 1)
    process.setArg("-nb", "gpu")
    process.start()
    process.wait()    
    print(process.getArgs())
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")
    return process


minimisation_steps = 500
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = int(sys.argv[1])

os.system(f"mkdir -p ../equilibration_files/ligands/")
os.system(f"mkdir -p ../equilibration_files/protein/")

ligand_stream = open("../ligands.dat", "r")
ligand_lines = ligand_stream.readlines()
ligand_name = ligand_lines[index].rstrip()

ligand_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_ligand_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_ligand_solvated.rst7"])
system_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_system_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_system_solvated.rst7"])

print("---------------------------")
print("Working on solvated ligand.")

print(f"\nMinimising in {minimisation_steps} steps.")
protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
minimised_process = run_process(ligand_solvated, protocol)

output_files = minimised_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/ligands/minimisation_{index}")

minimised = minimised_process.getSystem()

print("Finished minimisation.")


print(f"\nNVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = bss.Protocol.Equilibration(
                                    runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                    temperature_start=0*bss.Units.Temperature.kelvin,
                                    temperature_end=300*bss.Units.Temperature.kelvin,
                                    restraint="all"
                                )
restrained_nvt_process = run_process(minimised, protocol)

output_files = restrained_nvt_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/ligands/restrained_nvt_{index}")

restrained_nvt = restrained_nvt_process.getSystem()

print("Finished restrained NVT.")

print(f"\nNVT equilibration for {runtime_nvt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
nvt_process = run_process(restrained_nvt, protocol)

output_files = nvt_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/ligands/nvt_{index}")

nvt = restrained_nvt_process.getSystem()

print("Finished unrestrained NVT.")


print(f"\nNPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )
restrained_npt_process = run_process(nvt, protocol)

output_files = restrained_npt_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/ligands/restrained_npt_{index}")

restrained_npt = restrained_npt_process.getSystem()

print("Finished restrained NPT.")

print(f"\nNPT equilibration for {runtime_npt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )

npt_process = run_process(restrained_npt, protocol)

output_files = npt_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/ligands/npt_{index}")

equilibrated_ligand = npt_process.getSystem()

print("Finished equilibrating ligand.")

print("---------------------------")
print("Working on solvated ligand+protein.")

print(f"\nMinimising in {minimisation_steps} steps..")
protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
minimised_system_process = run_process(system_solvated, protocol)

output_files = minimised_system_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/minimisation_{index}")

minimised_system = minimised_system_process.getSystem()

print("Finished minimisation.")

print(f"\nNVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                temperature_start=0*bss.Units.Temperature.kelvin,
                                temperature_end=300*bss.Units.Temperature.kelvin,
                                restraint="all"
                                )
restrained_nvt_system_process = run_process(minimised_system, protocol)

output_files = restrained_nvt_system_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/restrained_nvt_{index}")

restrained_nvt_system = restrained_nvt_system_process.getSystem()

print("Finished restrained NVT.")

print(f"\nNVT equilibration for {runtime_nvt} ps while restraining all backbone atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="backbone"
                                )
backbone_restrained_nvt_process = run_process(restrained_nvt_system, protocol)

output_files = backbone_restrained_nvt_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/backbone_restrained_nvt_{index}")

backbone_restrained_nvt = backbone_restrained_nvt_process.getSystem()

print("Finished backbone-restrained NVT.")

print(f"\nNVT equilibration for {runtime_nvt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature_end=300*bss.Units.Temperature.kelvin,
                                )
nvt_system_process = run_process(backbone_restrained_nvt, protocol)

output_files = nvt_system_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/nvt_{index}")

nvt_system = nvt_system_process.getSystem()

print("Finished NVT.")

print(f"\nNPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )
restrained_npt_system_process = run_process(nvt_system, protocol)

output_files = restrained_npt_system_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/restrained_npt_{index}")

restrained_npt_system = restrained_npt_system_process.getSystem()

print("Finished restrained NPT.")

print(f"\nNPT equilibration for {runtime_npt} ps without restraints.")
protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
npt_system_process = run_process(restrained_npt_system, protocol)

output_files = npt_system_process.workDir()
os.system(f"cp -r {output_files} ../equilibration_files/protein/npt_{index}")

equilibrated_system = npt_system_process.getSystem()

print("Finished equilibrating ligand-protein system.")
    
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
