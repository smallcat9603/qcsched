#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=regular-o
#PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:02:00
#PJM -g gz00
#PJM -j
#------- Program execution -------#

module purge
module load fj
module load fjmpi
module load waitio
export WAITIO_MASTER_HOST=`hostname`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=1

mpiexec  ./test-a64fx
exit
