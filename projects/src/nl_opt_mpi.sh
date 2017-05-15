#!/bin/sh
#
# SGE batch script
# request Bourne shell as shell for job
#$ -N nl_opt_mpi
#$ -S /bin/bash
#$ -m e
##$ -M bengtsson@bnl.gov
#$ -cwd
#$ -o nl_opt_mpi_log.$JOB_ID -j y


echo "$JOB_ID.$TASK_ID starts at $HOSTNAME `date`"

dir=`pwd`
export SGE_O_WORKDIR=$dir

echo "Got $NSLOTS slots."
#cat $TMPDIR/machines

echo "# Working dir: $SGE_O_WORKDIR"
cd $SGE_O_WORKDIR
#ls -la

date

#unset SGE_ROOT
#unset LD_PRELOAD=/usr/local/bin/libmpich.so
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
mpiexec $HOME/Thor-2.0/thor/wrk/nl_opt_mpi 4

date

#qsub -q Parallel.q -l MPI=true -pe mpich2 <no of CPUs> nl_opt_mpi.sh

