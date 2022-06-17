#!/bin/bash

# General skeleton for running the FEP

# 1. Set the lambda values
# 2. Loop over transformation folders
# 3. Loop over bound and unbound folders
# 4. Loop over lambda folders
# 5. Do production run in each lambda folder using SOMD

lig0=$1
lig1=$2
GPU=$3

echo "OPENCL SOMD SCRIPT"

# Select which GPU to run on
export OPENMM_PLUGIN_DIR="/home/jguven/Software/miniconda3/envs/bss-d/lib/plugins/"
export OPENCL_VISIBLE_DEVICES=$GPU
export OPENMM_DEFAULT_PLATFORM="OpenCL"
echo "Running on device: $OPENCL_VISIBLE_DEVICES"

# Uncomment for testing
#lambdas=0.0000

lambdas=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )
engine="SOMD"

#IFS=',' read -r -a lambdas <<< "$lambdastring"
N_WINDOWS=${#lambdas[@]}

# Change directory to the perturbation folder
cd ../outputs/SOMD/$lig0~$lig1/
echo "cd ../outputs/SOMD/$lig0~$lig1/"

# Change into stage folder
cd bound
echo "cd bound"

# Loop over all the lambdas
for lambda in "${lambdas[@]}" 
do
echo "lambda is: " $lambda "at bound for " $lig0~$lig1

cd lambda_$lambda
echo "cd lambda_$lambda"
echo "production run..."

somd-freenrg -C ./somd.cfg -l $lambda -d $GPU -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p OpenCL 1> somd.log 2> somd.err
echo "somd-freenrg -C ./somd.cfg -l $lambda -d $GPU -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p OpenCL 1> somd.log 2> somd.err"
cd ../

echo "$PWD"
done

exit 0

