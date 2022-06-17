#!/bin/bash

N_LIGANDS=$(ls -l parallel_solvation/ligands | grep ^d | wc -l)

N_TASKS=$(($N_LIGANDS -1))


for i in $(seq 0 $N_TASKS)
do
#DIR_INDEX=$(($i +1))
echo "task index is $i"
python amber_to_gromacs.py $i
done

