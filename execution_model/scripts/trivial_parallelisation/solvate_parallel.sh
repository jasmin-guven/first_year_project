#!/bin/bash

N_LIGANDS=$(ls -l ligands | grep ^d | wc -l)

N_TASKS=$(($N_LIGANDS -1))

for i in $(seq 0 $N_TASKS)
do
DIR_INDEX=$(($i +1))
cd ligands/ligand_$DIR_INDEX
#nohup python solvation_$i.py $i > solv.out &
cd ../../

done
