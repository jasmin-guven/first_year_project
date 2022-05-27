import BioSimSpace as bss
import sys
from BioSimSpace import _Exceptions
import os


def run_process(system, md_protocol):
    """
    @param system: solvated system (BSS object)
    @param md_protocol: BSS protocolBSS protocol
    @return the processed system.
    """
    process = BSS.Process.Gromacs(system, md_protocol)
    process.setArg("-nt", 1)
    print(process.getOutput)
    process.start()
    process.wait()
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise _Exceptions.ThirdPartyError("The process exited with an error!")
    system = process.getSystem()
    return system


minimisation_steps = 250
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

print(f"program: {sys.argv[0]}, index: {sys.argv[1]}")
index = int(sys.argv[1])

ligand_stream = open("../../../../ligands.dat", "r")
ligand_lines = ligand_stream.readlines()
ligand_name = ligand_lines[index].rstrip()

ligand_solvated = BSS.IO.readMolecules([f"../inputs/ligands/{ligand_name}_ligand_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_ligand_solvated.rst7"])
system_solvated = BSS.IO.readMolecules([f"../inputs/ligands/{ligand_name}_system_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_system_solvated.rst7"])
print("---------------------------")
print("Working on solvated ligand.")
print(f"Minimising in {minimisation_steps} steps..")
protocol = BSS.Protocol.Minimisation(steps=minimisation_steps)
minimised = run_process(ligand_solvated, protocol)
print("Finished minimisation.")

print(f"NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = BSS.Protocol.Equilibration(
                                    runtime=runtime_short_nvt*BSS.Units.Time.picosecond,
                                    temperature_start=0*BSS.Units.Temperature.kelvin,
                                    temperature_end=300*BSS.Units.Temperature.kelvin,
                                    restraint="all"
                                )
restrained_nvt = run_process(minimised, protocol)

print(f"NVT equilibration for {runtime_nvt} ps without restraints.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_nvt*BSS.Units.Time.picosecond,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                )

nvt = run_process(restrained_nvt, protocol)

print(f"NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_npt*BSS.Units.Time.picosecond,
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                restraint="heavy",
                                )
restrained_npt = run_process(nvt, protocol)

print(f"NPT equilibration for {runtime_npt} ps without restraints.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_npt*BSS.Units.Time.picosecond,
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                )
equilibrated_ligand = run_process(restrained_npt, protocol)

print("---------------------------")
print("Working on solvated ligand+protein.")
print(f"Minimising in {minimisation_steps} steps..")
protocol = BSS.Protocol.Minimisation(steps=minimisation_steps)
minimised = run_process(system_solvated, protocol)

print(f"NVT equilibration for {runtime_short_nvt} ps while restraining all non-solvent atoms.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_short_nvt*BSS.Units.Time.picosecond,
                                temperature_start=0*BSS.Units.Temperature.kelvin,
                                temperature_end=300*BSS.Units.Temperature.kelvin,
                                restraint="all"
                                )
restrained_nvt_system = run_process(minimised, protocol)

print(f"NVT equilibration for {runtime_nvt} ps while restraining all backbone atoms.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_nvt*BSS.Units.Time.picosecond,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                restraint="backbone"
                                )
backbone_restrained_nvt_system = run_process(restrained_nvt_system, protocol)

print(f"NVT equilibration for {runtime_nvt} ps without restraints.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_nvt*BSS.Units.Time.picosecond,
                                temperature_end=300*BSS.Units.Temperature.kelvin,
                                )

nvt_system = run_process(backbone_restrained_nvt_system, protocol)

print(f"NPT equilibration for {runtime_npt} ps while restraining non-solvent heavy atoms..")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_npt*BSS.Units.Time.picosecond,
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                restraint="heavy",
                                )
restrained_npt_system = run_process(npt, protocol)

print(f"NPT equilibration for {runtime_npt} ps without restraints.")
protocol = BSS.Protocol.Equilibration(
                                runtime=runtime_npt*BSS.Units.Time.picosecond,
                                pressure=1*BSS.Units.Pressure.atm,
                                temperature=300*BSS.Units.Temperature.kelvin,
                                )
equilibrated_system = run_process(restrained_npt_system, protocol)

os.system("mkdir -p ../prep/ligands")
os.system("mkdir -p ../prep/protein")

print("Saving solvated/equilibrated systems.")
print("Ligand:")
print(equilibrated_ligand)
BSS.IO.saveMolecules(f"../prep/ligands/{ligand_name}_ligand_prepped", equilibrated_ligand, ["PRM7", "RST7"])

print("\n Ligand + protein:")
print(equilibrated_system)
BSS.IO.saveMolecules(f"../prep/protein/{lig_name}_system_prepped", equilibrated_system, ["PRM7", "RST7"])
print("First 20 molecules in ligand + protein system:")
for molecules in equilibrated_system.getMolecules()[:20]:
    print(molecules)
print("Done.")
