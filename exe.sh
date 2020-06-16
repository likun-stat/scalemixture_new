
#!/bin/bash

#PBS -A drh20_a_g_sc_default 
#PBS -V
#PBS -l nodes=1:ppn=20
#PBS -l walltime=48:00:00
#PBS -N scalemixture40
#PBS -j oe
#PBS -l pmem=8gb

#-#-# Echo
echo "#-#-#Job started on `hostname` at `date` "
echo This job runs on the following processors:
echo `cat $PBS_NODEFILE`


cd $PBS_O_WORKDIR
echo " "
echo " "
echo "Job started on hostname at date"


module purge
module load r/3.4

Rscript imputation_animation_02.R

echo "#-#-#Job Ended at `date`"
