#!/bin/sh
#PJM -N "s1"
#PJM -L rscgrp=lecture-a
#-PJM -L node=1
#PJM -L gpu=1
#PJM --mpi proc=1
#PJM --omp thread=1
#PJM -L elapse=00:15:00
#PJM -g gt00
#PJM -j
#PJM -e err
#PJM -o t1.lst

module unload aquarius
module unload gcc
module load intel
module load impi
module load metis

mpiexec -n ${PJM_MPI_PROC} ./sol1
