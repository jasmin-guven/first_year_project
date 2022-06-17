#!/bin/bash

mapfile PERTURBATIONS < ../network.dat
for perturbation in "${PERTURBATIONS[@]}"
do

IFS=' ' read -a line <<< "$perturbation"
ligand_1=${line[0]}
ligand_2=${line[1]}
n_windows=${line[2]}
lambdastring=${line[3]}
engine=${line[4]}

echo "@@@ Processing $ligand_1 $ligand_2 in $n_windows lambda windows @@@" 
nohup ./rbfe.sh $ligand_1 $ligand_2 $lambdastring $engine > nohup_output/somd_$ligand_2.out &
done
