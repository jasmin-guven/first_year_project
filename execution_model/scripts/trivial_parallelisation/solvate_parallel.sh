#!/bin/bash

DIRS=$(ls -l ligands | grep ^d | wc -l)

for i in $(seq 1 $DIRS)
do
cd ligands/ligand_$i
nohup python solvation_$i.py > solv.out &
cd ../../
done
