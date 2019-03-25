#!/bin/sh

prm1=${1-0}

queue="ap-high.q"
#queue="ap-medium.q"
#queue="test-medium.q"

#t1="00:10:00"
#t2="12:00:00"

dir=$THOR_LIB/projects/src

\rm tbs.cmd.o*

#qsub -l h_rt=$t2 -q $queue $dir/tbs.cmd
# NO = 9: ~30GB.
qsub -l mem_free=30G,h_vmem=30G -q $queue $dir/tbs.cmd
# NO = 11: ~50GB.
#qsub -l mem_free=60G,h_vmem=60G -q $queue $dir/tbs.cmd
#qsub -q $queue $dir/tbs.cmd
