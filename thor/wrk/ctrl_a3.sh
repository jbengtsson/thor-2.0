#!/bin/sh

dir=`pwd`

inf=$1

\rm ctrl_a3.cmd.o*

qsub -l Fast ctrl_a3.cmd $inf
