#!/bin/bash
cd equilibration/

ligands=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 )
# ligands=( 16 )
cd ligands
for ligand in "${ligands[@]}" 
do
cd lig_h_$ligand
echo $PWD
gmx energy -f minim.edr -o potential.xvg<<EOF
Potential
0
EOF
gmx energy -f nvt.edr -o temperature.xvg<<EOF
Temperature
0
EOF
gmx energy -f npt.edr -o pressure.xvg<<EOF
Pressure
0
EOF
gmx energy -f npt.edr -o density.xvg<<EOF
Density
0
EOF
cd ..
done
cd ../systems
for ligand in "${ligands[@]}" 
do
cd lig_h_$ligand
echo $PWD
gmx energy -f minim.edr -o potential.xvg<<EOF
Potential
0
EOF
gmx energy -f nvt.edr -o temperature.xvg<<EOF
Temperature
0
EOF
gmx energy -f npt.edr -o pressure.xvg<<EOF
Pressure
0
EOF
gmx energy -f npt.edr -o density.xvg<<EOF
Density
0
EOF
cd ..
done