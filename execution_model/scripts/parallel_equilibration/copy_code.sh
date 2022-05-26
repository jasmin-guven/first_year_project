#!/bin/bash

# get number of ligand directories
DIRS=$(ls -l ../parallel_solvation/ligands | grep ^d | wc -l)
#echo $DIRS

cd ligands

# loop over directories
for i in $(seq 1 $DIRS) 
do 
cd ligand_$i
cp ../../equilibration_parallel.py equilibration_$i.py
cd ..
done
