#!/bin/bash

reps=$1    #Specify reps as (no. of replicas - 1) eg: For 8 replicas, reps = 7

for i in $(seq 0 $reps)
do
    echo stat.$i
    mkdir plots_stat.$i
    cp visualize_restraint_scores.py plots_stat.$i/
    cd plots_stat.$i/
    for field in ConnectivityRestraint_None ExcludedVolumeSphere_None GaussianEMRestraint_None GaussianEMRestraint_None_CCC GaussianEMRestraint_sigma_None MinimumPairDistanceBindingRestraint0_Score_BLD10-SAS4 MinimumPairDistanceBindingRestraint10_Score_SAS4-SAS4 MinimumPairDistanceBindingRestraint1_Score_BLD10-SAS4 MinimumPairDistanceBindingRestraint2_Score_BLD10-SAS4 MinimumPairDistanceBindingRestraint3_Score_BLD10-SAS4 MinimumPairDistanceBindingRestraint4_Score_BLD10-BLD10 MinimumPairDistanceBindingRestraint5_Score_BLD10-BLD10 MinimumPairDistanceBindingRestraint6_Score_BLD10-BLD10 MinimumPairDistanceBindingRestraint7_Score_BLD10-BLD10 MinimumPairDistanceBindingRestraint8_Score_BLD10-BLD10 MinimumPairDistanceBindingRestraint9_Score_SAS4-SAS4 SinglePointMinGaussianRestraint0_Score_SPMGRBld10-0 SinglePointMinGaussianRestraint1_Score_SPMGRSas4-0 SinglePointMinGaussianRestraint2_Score_SPMGRBld10-1 SinglePointMinGaussianRestraint3_Score_SPMGRSas4-1
    do 
        echo $i
        echo $field
        $IMP python ~/imp-clean/imp/modules/pmi/pyext/src/process_output.py -f ../../run_1/stat.$i.out -t $field >> score_$field.txt
        python visualize_restraint_scores.py ./score_$field.txt $field
        rm score_$field.txt
    done
    rm visualize_restraint_scores.py
    cd ../
done
