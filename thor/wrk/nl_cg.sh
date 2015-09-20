#!/bin/sh

\rm nl_cg.cmd.o*

qsub -l Fast nl_cg.cmd
