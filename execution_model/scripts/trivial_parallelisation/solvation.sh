#!/bin/bash

# Number of ligands
N_LIGANDS=$(< "../ligands.dat" wc -l)

# Number of tasks
N_TASKS=$(($N_LIGANDS -1))
#echo $N_TASKS

# loop over number of tasks
for i in $(seq 0 $N_TASKS)
do
python solvation.py $i
#echo $i
done
