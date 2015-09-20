#!/bin/sh

\rm nl_pwl.cmd.o*

qsub -l Fast nl_pwl.cmd
