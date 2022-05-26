#!/bin/bash
cd ligands
for i in {1..16}
do
cp ligand_$i/lig_h_$i.mol2 .
done
