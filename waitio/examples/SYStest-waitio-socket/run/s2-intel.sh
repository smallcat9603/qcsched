#!/bin/sh
#PJM -N "SYStest-waitio-socket"
#PJM -L rscgrp=coupler-a
#PJM -L node=1
#PJM --mpi proc=32
#PJM -L elapse=00:15:00
#PJM -g gz00
#PJM -j
#PJM -e err

module purge
module load intel
module load impi
module load metis
module load waitio

export WAITIO_MASTER_HOST=`waitio-serv -c`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=1
export WAITIO_NPB=3

hostname


mkdir s2; cp s2.dat s2/
(cd s2; mpiexec -n ${PJM_MPI_PROC} ../sol2-intel; cp s2.out ../; rm work.* )
rm -rf s2
#rm work.* 
