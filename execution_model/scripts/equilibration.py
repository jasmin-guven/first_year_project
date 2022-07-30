import BioSimSpace as bss
import numpy as np
import os


def change_barostat(system: bss._SireWrappers._system.System, 
                    protocol: bss.Protocol._equilibration.Equilibration, 
                    work_directory: str, 
                    process_name: str) -> None:
    """
    Change barostat in .mdp file for NPT runs
    @param: system
    @param: protocol
    @param: work_directory
    @param: process name
    @return: None
    """
    try: 
        process = bss.Process.Gromacs(system,
                                      protocol,
                                      name=process_name,
                                      work_dir=work_directory)
    except RuntimeError:

        with open(f"{work_directory}/{process_name}.mdp", "r") as mdp:
            lines = mdp.read()

        new_lines = lines.replace("pcoupl = berendsen", "pcoupl = C-rescale")    
        with open(f"{work_directory}/{process_name}.mdp", "w") as mdp:
            mdp.write(new_lines)
                

minimisation_steps = 1000
runtime_short_nvt = 5  # ps
runtime_nvt = 50  # ps
runtime_npt = 200  # ps

indices = np.arange(0, 16, 1)
# indices = [0]
for index in indices:
    ligand_stream = open("../ligands.dat", "r")
    ligand_lines = ligand_stream.readlines()
    ligand_name = ligand_lines[index].rstrip()
    print(f"Working on {ligand_name}.")
    print("=======================================")
    print("==============FREE STAGE===============")
    print("=======================================")
    os.system(f"mkdir -p equilibration/ligands/{ligand_name}")
    lig_work_dir = f"equilibration/ligands/{ligand_name}"
    ligand_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_ligand_solvated.prm7",
                      f"../inputs/ligands/{ligand_name}_ligand_solvated.rst7"])
    print("---------------------------------------")
    print("-------------MINIMISATION--------------")
    print("---------------------------------------")
    minim_lig_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
    minimisation_process = bss.Process.Gromacs(ligand_solvated, 
                                               minim_lig_protocol, 
                                               name="minim",
                                               work_dir=lig_work_dir)
    os.system(f"gmx grompp -f {lig_work_dir}/minim.mdp -c {lig_work_dir}/minim.gro -p {lig_work_dir}/minim.top -o {lig_work_dir}/minim.tpr")
    os.system(f"gmx mdrun -v -deffnm {lig_work_dir}/minim -nt 1 -nb gpu")bi
    print("---------------------------------------")
    print("-------------RESTRAINED NVT------------")
    print("---------------------------------------")
    minim_ligand = bss.IO.readMolecules([f"{lig_work_dir}/minim.gro",
                                         f"{lig_work_dir}/minim.top"])
    rnvt_lig_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                          temperature_start=0*bss.Units.Temperature.kelvin,
                                          temperature_end=300*bss.Units.Temperature.kelvin,
                                          restraint="all")

    rnvt_lig_process = bss.Process.Gromacs(minim_ligand, 
                                                 rnvt_lig_protocol, 
                                                 name="r_nvt",
                                                 work_dir=lig_work_dir)
    os.system(f"gmx grompp -f {lig_work_dir}/r_nvt.mdp -c {lig_work_dir}/minim.gro -r {lig_work_dir}/minim.gro -p {lig_work_dir}/r_nvt.top -o {lig_work_dir}/r_nvt.tpr")
    os.system(f"gmx mdrun -v -deffnm {lig_work_dir}/r_nvt -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("------------UNRESTRAINED NVT-----------")
    print("---------------------------------------")
    rnvt_ligand = bss.IO.readMolecules([f"{lig_work_dir}/r_nvt.gro",
                                        f"{lig_work_dir}/r_nvt.top"])
    nvt_lig_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
    nvt_lig_process = bss.Process.Gromacs(rnvt_ligand, 
                                          nvt_lig_protocol, 
                                          name="nvt",
                                          work_dir=lig_work_dir)
    os.system(f"gmx grompp -f {lig_work_dir}/nvt.mdp -c {lig_work_dir}/r_nvt.gro -p {lig_work_dir}/nvt.top -t {lig_work_dir}/r_nvt.cpt -o {lig_work_dir}/nvt.tpr")
    os.system(f"gmx mdrun -v -deffnm {lig_work_dir}/nvt -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("-------------RESTRAINED NPT------------")
    print("---------------------------------------")
    nvt_ligand = bss.IO.readMolecules([f"{lig_work_dir}/nvt.gro",
                                       f"{lig_work_dir}/nvt.top"])
    rnpt_lig_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )

    change_barostat(nvt_ligand, rnpt_lig_protocol, lig_work_dir, "r_npt")

    os.system(f"gmx grompp -f {lig_work_dir}/r_npt.mdp -c {lig_work_dir}/nvt.gro -r {lig_work_dir}/nvt.gro -p {lig_work_dir}/r_npt.top -t {lig_work_dir}/nvt.cpt -o {lig_work_dir}/r_npt.tpr")
    os.system(f"gmx mdrun -v -deffnm {lig_work_dir}/r_npt -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("------------UNRESTRAINED NPT-----------")
    print("---------------------------------------")
    npt_ligand = bss.IO.readMolecules([f"{lig_work_dir}/r_npt.gro",
                                       f"{lig_work_dir}/r_npt.top"])
    npt_lig_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin
                                )
    change_barostat(npt_ligand, npt_lig_protocol, lig_work_dir, "npt")
    os.system(f"gmx grompp -f {lig_work_dir}/npt.mdp -c {lig_work_dir}/r_npt.gro -t {lig_work_dir}/r_npt.cpt -p {lig_work_dir}/npt.top -o {lig_work_dir}/npt.tpr")
    os.system(f"gmx mdrun -v -deffnm {lig_work_dir}/npt -nt 1 -nb gpu")
    
    print("=======================================")
    print("==============BOUND STAGE==============")
    print("=======================================")
    os.system(f"mkdir -p equilibration/systems/{ligand_name}")
    sys_work_dir = f"equilibration/systems/{ligand_name}"
    system_solvated = bss.IO.readMolecules([f"../inputs/ligands/{ligand_name}_system_solvated.prm7",
                                        f"../inputs/ligands/{ligand_name}_system_solvated.rst7"])

    print("---------------------------------------")
    print("-------------MINIMISATION--------------")
    print("---------------------------------------")
    minim_sys_protocol = bss.Protocol.Minimisation(steps=minimisation_steps)
    minimisation_process = bss.Process.Gromacs(system_solvated, 
                                               minim_sys_protocol, 
                                               name="minim",
                                               work_dir=sys_work_dir)
    os.system(f"gmx grompp -f {sys_work_dir}/minim.mdp -c {sys_work_dir}/minim.gro -p {sys_work_dir}/minim.top -o {sys_work_dir}/minim.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/minim -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("-------------RESTRAINED NVT------------")
    print("---------------------------------------")
    minim_system = bss.IO.readMolecules([f"{sys_work_dir}/minim.gro",
                                         f"{sys_work_dir}/minim.top"])
    rnvt_sys_protocol = bss.Protocol.Equilibration(runtime=runtime_short_nvt*bss.Units.Time.picosecond,
                                          temperature_start=0*bss.Units.Temperature.kelvin,
                                          temperature_end=300*bss.Units.Temperature.kelvin,
                                          restraint="all")

    rnvt_sys_process = bss.Process.Gromacs(minim_system, 
                                           rnvt_sys_protocol, 
                                           name="r_nvt",
                                           work_dir=sys_work_dir)
    os.system(f"gmx grompp -f {sys_work_dir}/r_nvt.mdp -c {sys_work_dir}/minim.gro -r {sys_work_dir}/minim.gro -p {sys_work_dir}/r_nvt.top -o {sys_work_dir}/r_nvt.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/r_nvt -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("-----------BB RESTRAINED NVT-----------")
    print("---------------------------------------")
    rnvt_system = bss.IO.readMolecules([f"{sys_work_dir}/r_nvt.gro",
                                        f"{sys_work_dir}/r_nvt.top"])
    bbnvt_sys_protocol = bss.Protocol.Equilibration(runtime=runtime_nvt*bss.Units.Time.picosecond,
                                                    temperature_end=300*bss.Units.Temperature.kelvin,
                                                    restraint="backbone")
    bbnvt_sys_process = bss.Process.Gromacs(rnvt_system,
                                            bbnvt_sys_protocol,
                                            name="bb_nvt",
                                            work_dir=sys_work_dir)
    os.system(f"gmx grompp -f {sys_work_dir}/bb_nvt.mdp -c {sys_work_dir}/r_nvt.gro -r {sys_work_dir}/r_nvt.gro -p {sys_work_dir}/bb_nvt.top -t {sys_work_dir}/r_nvt.cpt -o {sys_work_dir}/bb_nvt.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/bb_nvt -nt 1 -nb gpu")

    print("---------------------------------------")
    print("------------UNRESTRAINED NVT-----------")
    print("---------------------------------------")
    bbnvt_system = bss.IO.readMolecules([f"{sys_work_dir}/bb_nvt.gro",
                                         f"{sys_work_dir}/bb_nvt.top"])
    nvt_sys_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_nvt*bss.Units.Time.picosecond,
                                temperature=300*bss.Units.Temperature.kelvin,
                                )
    nvt_sys_process = bss.Process.Gromacs(bbnvt_system, 
                                          nvt_sys_protocol, 
                                          name="nvt",
                                          work_dir=sys_work_dir)
    os.system(f"gmx grompp -f {sys_work_dir}/nvt.mdp -c {sys_work_dir}/bb_nvt.gro -p {sys_work_dir}/nvt.top -t {sys_work_dir}/bb_nvt.cpt -o {sys_work_dir}/nvt.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/nvt -nt 1 -nb gpu")
    
    
    print("---------------------------------------")
    print("-------------RESTRAINED NPT------------")
    print("---------------------------------------")
    nvt_system = bss.IO.readMolecules([f"{sys_work_dir}/nvt.gro",
                                       f"{sys_work_dir}/nvt.top"])
    rnpt_sys_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin,
                                restraint="heavy",
                                )

    change_barostat(nvt_system, rnpt_sys_protocol, sys_work_dir, "r_npt")

    os.system(f"gmx grompp -f {sys_work_dir}/r_npt.mdp -c {sys_work_dir}/nvt.gro -r {sys_work_dir}/nvt.gro -p {sys_work_dir}/r_npt.top -t {sys_work_dir}/nvt.cpt -o {sys_work_dir}/r_npt.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/r_npt -nt 1 -nb gpu")
    
    print("---------------------------------------")
    print("------------UNRESTRAINED NPT-----------")
    print("---------------------------------------")
    npt_system = bss.IO.readMolecules([f"{sys_work_dir}/r_npt.gro",
                                       f"{sys_work_dir}/r_npt.top"])
    npt_sys_protocol = bss.Protocol.Equilibration(
                                runtime=runtime_npt*bss.Units.Time.picosecond,
                                pressure=1*bss.Units.Pressure.atm,
                                temperature=300*bss.Units.Temperature.kelvin
                                )
    change_barostat(npt_system, npt_sys_protocol, sys_work_dir, "npt")
    os.system(f"gmx grompp -f {sys_work_dir}/npt.mdp -c {sys_work_dir}/npt.gro -t {sys_work_dir}/r_npt.cpt -p {sys_work_dir}/npt.top -o {sys_work_dir}/npt.tpr")
    os.system(f"gmx mdrun -v -deffnm {sys_work_dir}/npt -nt 1 -nb gpu")
  

    print("=======================================")
    print("=================DONE==================")
    print("=======================================")
