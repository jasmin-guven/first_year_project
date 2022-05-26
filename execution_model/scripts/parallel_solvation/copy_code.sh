#!/bin/bash

# get number of ligand directories
DIRS=$(ls -l ligands | grep ^d | wc -l)
#echo $DIRS

cd ligands

# loop over directories
for i in $(seq 1 $DIRS) 
do 
cd ligand_$i
cp ../../solvation_parallel.py solvation_$i.py
cd ..
done
