#!/bin/sh

\rm nl_ns.cmd.o*

qsub -l Fast nl_ns.cmd
