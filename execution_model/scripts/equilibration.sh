#!/bin/bash

# Number of ligands
n_ligands=$(< "../ligands.dat" wc -l)

# Number of tasks
n_tasks=$(($n_ligands -1))

# loop over number of tasks
for i in n_tasks:
do
	python equilibration.py $i
done
