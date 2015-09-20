#!/bin/sh

dir=`pwd`

~/projects/src/main /home/bengtsson/projects/in/lattice/mv15-c5-v2-May14-tracy 2
wait

\rm get_dnu.cmd.o* fmap*.cmd.o*

qsub -l Fast get_dnu.cmd
qsub -l Fast fmap.cmd
qsub -l Fast fmap_dp.cmd
