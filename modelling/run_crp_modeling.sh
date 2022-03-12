#!/bin/bash

# Usage:
# ./run_rnapolii_modeling.py out_dir n n_steps

py=/bin/python
NRUNS=$1    #no. of replicas to start with

#/home/kartik/imp-clean/build/setup_environment.sh 
for runid in $(seq $NRUNS)  
do 
    mpirun -np 8 $IMP python crp_modeling_final.py $runid &
done


