#!/bin/sh

prm1=${1-0}

queue="ap-high.q"
#queue="ap-medium.q"
#queue="test-medium.q"

#t1="00:10:00"
#t2="12:00:00"

dir=$THOR_LIB/projects/src

cd $prm1/dnu

\rm tune_conf_2.cmd.o*

#qsub -l h_rt=$t2 -q $queue $dir/sls_2.cmd
# NO = 11: ~50GB.
#qsub -l mem_free=60G,h_vmem=60G -q $queue -v ns=$prm1 $dir/tune_conf_2.cmd
qsub -q $queue -v ns=$prm1 $dir/tune_conf_2.cmd
