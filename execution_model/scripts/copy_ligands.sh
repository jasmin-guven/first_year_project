#!/bin/bash
cd ../../kpc2/ligands
for i in {1..16}
do
cp ligand_$i/lig_h_$i.mol2 .
done
