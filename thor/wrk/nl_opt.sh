#!/bin/sh

\rm nl_opt_mpi_log.*

qsub -q Parallel.q -l MPI=true -pe mpich2 10 nl_opt_mpi.sh
