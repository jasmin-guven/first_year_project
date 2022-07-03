c#!/bin/bash

N_LIGANDS=$(ls -l parallel_solvation/ligands | grep ^d | wc -l)

for i in $(seq 1 $N_LIGANDS)
do
#echo $i
PERTURBATION_INDEX=$(($i +1))
if [ $PERTURBATION_INDEX -lt 17 ]
then
echo "perturbation index is $PERTURBATION_INDEX"
nohup python rbfe_setup.py 1 $PERTURBATION_INDEX > nohup_output/setup_$i.out &
fi
done
