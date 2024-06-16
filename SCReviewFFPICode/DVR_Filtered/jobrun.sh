#PBS -N MQCIVR
#PBS -l nodes=1:ppn=12,walltime=10:00:00
#PBS -S /bin/bash
#PBS -q test
echo $PBS_NODEFILE
echo `cat $PBS_NODEFILE`

export GMPICONF=nodeinfo/PBS_JOBID

cd $PBS_O_WORKDIR
cp * $TMPDIR/
cd $TMPDIR
time ./dvr.x
mpirun -np 12 ./dyn.x
cp -f $TMPDIR/* $PBS_O_WORKDIR/


