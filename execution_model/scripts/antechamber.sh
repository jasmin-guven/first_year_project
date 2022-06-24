#!/bin/bash
cd ../../kpc2/ligands/
for i in {1..16}
do
mkdir ligand_$i
echo $PWD
mv lig_h_$i.mol2 ligand_$i/
cd ligand_$i/
echo "antechamber $PWD"
antechamber -at 2 -c bcc -s 2 -nc -1 -i lig_$i.mol2 -fi mol2 -o lig_h_$i.mol2 -fo mol2
cd ..
done
