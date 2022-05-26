#!/bin/bash
N_LIGANDS=$(< "../ligands.dat" wc -l)
DIRS=$(<"ligands/" grep -c ^d)
echo $DIRS

#cd ligands
#DIRS=$(($ligands))
#echo $DIRS
#for i in {1..15} #loop over 1-15 
#do 
#cd id_file_$i #cd into the new directories
#cp ../../retrieve_pdbs_parallel.py retrieve_pdbs_parallel_$i.py #copy code into each folder
#cd ../
#done
