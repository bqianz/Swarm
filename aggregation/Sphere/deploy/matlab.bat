#!/bin/bash
#PBS -S /bin/bash

# Choose the MCR directory according to the compiler version used
MCR=/global/software/matlab/mcr/v90

echo "Running on host: `hostname`"
cd $PBS_O_WORKDIR 
echo "Current working directory is `pwd`"

echo "Starting run at: `date`" 
./main_spherical.exe $MCR > mycode_${PBS_JOBID}.out
echo "Job finished at: `date`"