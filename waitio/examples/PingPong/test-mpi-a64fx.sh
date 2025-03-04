#!/bin/sh
#------ pjsub option --------#
#PJM -N "test1-mpi"
#PJM -L rscgrp=lecture-o
#PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:02:00
#PJM -g gt00
#PJM -j
#------- Program execution -------#

hostname
module purge
module load fj
module load fjmpi
module load waitio
export WAITIO_MASTER_HOST=`hostname`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=1

#waitio-serv-a64fx -m $WAITIO_MASTER_HOST
echo  "serv host set ${WAITIO_MASTER_HOST} "

mpiexec  ./test-mpi-a64fx
exit
