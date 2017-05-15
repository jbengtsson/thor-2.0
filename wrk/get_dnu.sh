#!/bin/sh

dir=`pwd`

~/projects/src/main commissioning_symm 5
wait

\rm get_dnu.cmd.o* fmap*.cmd.o*

qsub -l Fast get_dnu.cmd
qsub -l Fast fmap.cmd
qsub -l Fast fmap_dp.cmd
