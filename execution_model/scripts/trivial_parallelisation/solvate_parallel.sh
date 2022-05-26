#!/bin/bash

N_LIGANDS=$(ls -l ligands | grep ^d | wc -l)

N_TASKS=$(($N_LIGANDS -1))

for i in $(seq 0 $N_TASKS)
do

cd ligands/ligand_$i
nohup python solvation_$i.py $i > solv.out &
cd ../../

done
