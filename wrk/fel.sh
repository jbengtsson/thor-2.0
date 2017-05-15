#!/bin/sh

dir=`pwd`

\rm fel.cmd.o*

qsub -l Fast fel.cmd
