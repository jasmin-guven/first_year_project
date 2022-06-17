#!/bin/bash

# General skeleton for running the FEP

# 1. Set the lambda values
# 2. Loop over transformation folders
# 3. Loop over bound and unbound folders
# 4. Loop over lambda folders
# 5. Do production run in each lambda folder using SOMD

echo "OPENCL SOMD SCRIPT"

# Select which GPU to run on
export OPENMM_PLUGIN_DIR="/home/jguven/Software/miniconda3/envs/bss-d/lib/plugins/"
export OPENCL_VISIBLE_DEVICES=0
export OPENMM_DEFAULT_PLATFORM="OpenCL"
echo "Running on device: $OPENCL_VISIBLE_DEVICES"

lig0=$1
lig1=$2

# Uncomment for testing
lambdas=(0.0000 0.1000)

#lambdas=( 0.0000 0.1000 0.2000 0.3000 0.4000 0.5000 0.6000 0.7000 0.8000 0.9000 1.0000 )
engine="SOMD"

#IFS=',' read -r -a lambdas <<< "$lambdastring"
N_WINDOWS=${#lambdas[@]}

# Change directory to the perturbation folder
cd ../outputs/SOMD/test/
echo "cd ../outputs/SOMD/test/"

for stage in bound 
do

# Change into stage folder
cd $stage
echo "cd $stage"

# Loop over all the lambdas
for lambda in "${lambdas[@]}" 
do
echo "lambda is: " $lambda "at" $stage " for " $lig0~$lig1

cd lambda_$lambda
echo "$cd lambda_$lambda"
echo "production run..."

somd-freenrg -C ./somd.cfg -l $lambda -c ./somd.rst7 -t ./somd.prm7 -m ./somd.pert -p OpenCL
cd ..
echo "cd .."
done
cd ..
echo "cd .."
done
exit 0
