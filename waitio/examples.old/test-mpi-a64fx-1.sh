#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=debug-o
#PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:02:00
#PJM -g jh210022o
#PJM -j
#------- Program execution -------#

hostname
export LD_LIBRARY_PATH=/work/share/waitio/lib/:$LD_LIBRARY_PATH
export PATH=/work/share/waitio/bin/:$PATH
module load fj
module load fjmpi
export WAITIO_MASTER_HOST=`hostname`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=2

mpiexec  ./test-mpi-a64fx
exit
