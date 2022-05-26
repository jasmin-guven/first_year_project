#!/bin/bash

# Number of ligands
n_ligands=$(< "../ligands.dat" wc -l)

# Number of tasks
n_tasks=$(( $n_ligands -1 ))

# Perturbations
mapfile PERTURBATIONS < ../network.dat
#echo ${PERTURBATIONS[@]}

for ligand


#for perturbation in "${PERTURBATIONS[@]}"
#do
#	#echo $perturbation
	#IFS=' ' read -a ligand_pair <<< "$perturbation"
	#echo $ligand_pair
#done

