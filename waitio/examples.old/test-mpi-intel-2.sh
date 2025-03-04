#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=debug-a
#PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:05:00
#PJM -g jh210022a
#PJM -j
#------- Program execution -------#

module unload aquarius
module unload gcc
module load intel
module load impi

export LD_LIBRARY_PATH=/work/share/waitio/lib/:$LD_LIBRARY_PATH
export PATH=/work/share/waitio/bin/:$PATH

#export WAITIO_MASTER_HOST=wo0001-ib0
export WAITIO_MASTER_HOST=wo0009
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=1
export WAITIO_NPB=2

mpiexec -n 2  ./test-mpi

exit
