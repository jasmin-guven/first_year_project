#!/bin/bash

mkdir lig_minimisation
mkdir lig_npt
mkdir lig_nvt
mkdir sys_minimisation
mkdir sys_npt
mkdir sys_nvt

mapfile PERTURBATIONS < ../network_bound.dat
for perturbation in "${PERTURBATIONS[@]}"
do

#echo $perturbation

INPUTFILESTREAM=' ' read -a line <<< "$perturbation"

ligand_2=${line[1]}
echo $ligand_2
ligand_gro=../inputs/ligands/${ligand_2}_ligand_solvated.gro
ligand_top=../inputs/ligands/${ligand_2}_ligand_solvated.top
system_gro=../inputs/ligands/${ligand_2}_system_solvated.gro
system_top=../inputs/ligands/${ligand_2}_system_solvated.top

echo $ligand_gro $system_gro

mkdir lig_minimisation/lig_h_1~${ligand_2}
mkdir lig_npt/lig_h_1~${ligand_2}
mkdir lig_nvt/lig_h_1~${ligand_2}
mkdir sys_minimisation/lig_h_1~${ligand_2}
mkdir sys_npt/lig_h_1~${ligand_2}
mkdir sys_nvt/lig_h_1~${ligand_2}

echo "Working on solvated ligand: $ligand_2"
echo "Minimising ligand"

gmx grompp -f ../inputs/gromacs_files/ligand/minimisation.mdp -c $ligand_gro -p $ligand_top -o lig_minimisation/lig_h_1~${ligand_2}/min.tpr 
gmx mdrun -deffnm lig_minimisation/lig_h_1~${ligand_2}/min -nt 1 -nb gpu

echo "NVT equilibration for 5 ps while restraining all non-solvent atoms"
gmx grompp -f ../inputs/gromacs_files/ligand/restrained_nvt.mdp -c $ligand_gro -p $ligand_top -o lig_nvt/lig_h_1~${ligand_2}/restrained_nvt.tpr
gmx mdrun -deffnm lig_nvt/lig_h_1~${ligand_2}/restrained_nvt -nt 1 -nb gpu 
echo "Finished restrained NVT"

echo "NVT equilibration for 50 ps without restraints"
gmx grompp -f ../inputs/gromacs_files/ligand/nvt.mdp -c $ligand_gro -p $ligand_top -o lig_nvt/lig_h_1~${ligand_2}/nvt.tpr
gmx mdrun -deffnm lig_nvt/lig_h_1~${ligand_2}/nvt -nt 1 -nb gpu
echo "Finished NVT"

echo "NPT equilibration for 200 ps while restraining non-solvent heavy atoms"
gmx grompp -f ../inputs/gromacs_files/ligand/restrained_npt.mdp -c $ligand_gro -p $ligand_top -o lig_npt/lig_h_1~${ligand_2}/restrained_npt.tpr
gmx mdrun -deffnm lig_npt/lig_h_1~${ligand_2}/restrained_npt -nt 1 -nb gpu
echo "Finished restrained NPT"

echo "NPT equilibration for 200 ps without restraints"
gmx grompp -f ../inputs/gromacs_files/ligand/npt.mdp -c $ligand_gro -p $ligand_top -o lig_npt/lig_h_1~${ligand_2}/npt.tpr
gmx mdrun -deffnm lig_npt/lig_h_1~${ligand_2}/npt -nt 1 -nb gpu
echo "Finished NPT"

echo "Working on solvated ligand + protein system: $ligand_2"
echo "Minimising system"

gmx grompp -f ../inputs/gromacs_files/system/minimisation.mdp -c $system_gro -p $system_top -o sys_minimisation/lig_h_1~${ligand_2}/min.tpr 
gmx mdrun -deffnm sys_minimisation/lig_h_1~${ligand_2}/min -nt 1 -nb gpu
c
echo "NVT equilibration for 5 ps while restraining all non-solvent atoms"
gmx grompp -f ../inputs/gromacs_files/system/restrained_nvt.mdp -c $system_gro -p $system_top -o sys_nvt/lig_h_1~${ligand_2}/restrained_nvt.tpr
gmx mdrun -deffnm sys_nvt/lig_h_1~${ligand_2}/restrained_nvt -nt 1 -nb gpu 
echo "Finished restrained NVT"

echo "NVT equilibration for 50 ps while restraining all backbone atoms"
gmx grompp -f ../inputs/gromacs_files/system/bb_restrained_nvt.mdp -c $system_gro -p $system_top -o sys_nvt/lig_h_1~${ligand_2}/bb_restrained_nvt.tpr
gmx mdrun -deffnm sys_nvt/lig_h_1~${ligand_2}/bb_restrained_nvt -nt 1 -nb gpu
echo "Finished backbone-restrained NVT"

echo "NVT equilibration for 50 ps without restraints"
gmx grompp -f ../inputs/gromacs_files/system/nvt.mdp -c $system_gro -p $system_top -o sys_nvt/lig_h_1~${ligand_2}/nvt.tpr
gmx mdrun -deffnm sys_nvt/lig_h_1~${ligand_2}/nvt -nt 1 -nb gpu
echo "Finished NVT"

echo "NPT equilibration for 200 ps while restraining non-solvent heavy atoms"
gmx grompp -f ../inputs/gromacs_files/system/restrained_npt.mdp -c $system_gro -p $system_top -o sys_npt/lig_h_1~${ligand_2}/restrained_npt.tpr
gmx mdrun -deffnm sys_npt/lig_h_1~${ligand_2}/restrained_npt -nt 1 -nb gpu
echo "Finished restrained NPT"

echo "NPT equilibration for 200 ps without restraints"
gmx grompp -f ../inputs/gromacs_files/system/npt.mdp -c $system_gro -p $system_top -o sys_npt/lig_h_1~${ligand_2}/npt.tpr
gmx mdrun -deffnm sys_npt/lig_h_1~${ligand_2}/npt -nt 1 -nb gpu
echo "Finished NPT"
mv \#*\# backups/
mv step* backups
done
