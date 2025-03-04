#!/bin/sh
#PJM -N "SYStest-waitio-socket"
#PJM -L rscgrp=coupler-a
#PJM -L node=1
#PJM --mpi proc=1
#PJM --omp thread=1
#PJM -L elapse=00:15:00
#PJM -g gz00
#PJM -j
#PJM -e err

module purge
module load intel
module load impi
module load metis
module load waitio

export WAITIO_MASTER_HOST=`hostname`-ib0
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=3

hostname
waitio-serv -d -m ${WAITIO_MASTER_HOST}

echo $WAITIO_MASTER_HOST
mpiexec -n ${PJM_MPI_PROC} ./sol1-intel
#rm work.* aaa0.*
