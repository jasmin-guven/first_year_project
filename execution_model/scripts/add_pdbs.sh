#!/bin/bash

ligs=( 2 3 4 5 6 7 8 9 10 11 13 14 15 16 )

for lig in "${ligs[@]}"
do
cp ../outputs/SOMD/lig_h_1~lig_h_$lig/bound/lambda_1.0000/latest.pdb ../outputs/SOMD/PDBs/latest_$lig.pdb
cp ../outputs/SOMD/lig_h_1~lig_h_$lig/bound/lambda_1.0000/somd.prm7 ../outputs/SOMD/PDBs/somd_$lig.prm7
done
 
