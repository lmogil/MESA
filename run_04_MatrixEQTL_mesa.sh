#!/bin/bash
#PBS -N testjob
#PBS -S /bin/bash
#PBS -l walltime=100:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=10gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd /home/lauren/scripts

Rscript 04_MatrixEQTL_mesa.R --args ALL 10 
Rscript 04_MatrixEQTL_mesa.R --args ALL 20 
Rscript 04_MatrixEQTL_mesa.R --args ALL 30 
Rscript 04_MatrixEQTL_mesa.R --args ALL 50 
Rscript 04_MatrixEQTL_mesa.R --args ALL 100 
