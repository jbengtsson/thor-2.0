#!/bin/sh

#queue="ap-high.q"
queue="ap-medium.q"
#queue="test-medium.q"

#t1="00:10:00"
#t2="12:00:00"

dir=$THOR_LIB/projects/src

\rm diamond_ii.cmd.o*

#qsub -l h_rt=$t2 -q $queue $dir/sls_2.cmd
qsub -q $queue $dir/diamond_ii.cmd
