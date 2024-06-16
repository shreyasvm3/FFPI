#!/bin/bash

#SBATCH -J gridYZ    
#SBATCH -t 168:00:00
#SBATCH -n 12
#SBATCH -N 1
#SBATCH -p astra
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH --exclusive 
#SBATCH -w c0055 
# Enter the working directory
cd ${SLURM_SUBMIT_DIR}

echo "starting $SLURM_JOB_ID at `date` on `hostname`"

echo "$USER"

/usr/bin/mkdir -p /tmp/$USER/$SLURM_JOB_ID
#echo "Copying data over... "
cp $SLURM_SUBMIT_DIR/dvr.x /tmp/$USER/$SLURM_JOB_ID/. 
cp $SLURM_SUBMIT_DIR/dyn.x /tmp/$USER/$SLURM_JOB_ID/. 
cp $SLURM_SUBMIT_DIR/input /tmp/$USER/$SLURM_JOB_ID/. 

cd /tmp/$USER/$SLURM_JOB_ID/

#run fortran code
time ./dvr.x
time mpiexec -np 12 ./dyn.x

# Copy output files back from the temp directory to working directory
rsync -r /tmp/$USER/$SLURM_JOB_ID/ $SLURM_SUBMIT_DIR/

rm -rf /tmp/$USER/$SLURM_JOB_ID/
exit 0 


 
