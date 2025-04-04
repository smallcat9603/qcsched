#!/bin/sh
# This script is for multiple PB script test for 2 PB (1/2)
#PJM -N "SYStest-waitio-socket"
#PJM -L rscgrp=coupler-o
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1
#PJM -L elapse=00:15:00
#PJM -g gz00
#PJM -j
#PJM -e err

module purge
module load fj
module load fjmpi
module load metis
module load waitio

export WAITIO_MASTER_HOST=`hostname`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=3

hostname
waitio-serv-a64fx -d -m $WAITIO_MASTER_HOST

mpiexec -n ${PJM_MPI_PROC} ./sol1-a64fx
#rm work.* aaa0.*
