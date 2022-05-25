#!/bin/bash
cd ligands
for i in {1..16}
do 
cd ligand_$i/
parmchk2 -i lig_h_$i.mol2 -o lig_$i.frcmod -f mol2 -s 2 -a Y -w Y
cd ..
done
