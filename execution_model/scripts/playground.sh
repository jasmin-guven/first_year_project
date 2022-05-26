#!/bin/bash

# Number of ligands
n_ligands=$(< "../ligands.dat" wc -l)

# Number of tasks
n_tasks=$(( $n_ligands -1 ))

# Perturbations
mapfile perturbations < ../network.dat
echo $perturbations

