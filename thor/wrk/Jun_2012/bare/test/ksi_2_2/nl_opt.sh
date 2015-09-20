#!/bin/sh

dir=`pwd`

\rm nl_opt.cmd.o*

qsub -l Fast nl_opt.cmd
