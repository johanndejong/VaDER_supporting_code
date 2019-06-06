#!/bin/bash

source ~/.bash_profile

DATE=`date +"%Y%m%d%H%M%S"`
tempdir=PPMI_hyperparameter_optimization
mkdir $tempdir

for i in {1..20}
do
	sleep 1 # make sure that CPU load is spread in time
	Rscript PPMI_hyperparameter_optimization.r 6 ${i} ${DATE} > ${tempdir}/${i}.out 2>&1 &
done

