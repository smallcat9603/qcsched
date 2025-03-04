#!/bin/sh
#PJM -N "s3"
#PJM -L rscgrp=debug-a
#PJM -L node=1
#PJM --mpi proc=32
#PJM -L elapse=00:30:00
#- PJM -g jh180023
#PJM -g jh210022a
#PJM -j
#PJM -e err
#PJM -o t3.lst

module unload aquarius
module unload gcc
module load intel
module load impi
module load metis

mpiexec -n ${PJM_MPI_PROC} ./sol3

rm wk*
