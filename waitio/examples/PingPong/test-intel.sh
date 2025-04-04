#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=lecture-a
#-PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:05:00
#PJM -L gpu=1
#PJM -g gt00
#PJM -j
#------- Program execution -------#

module purge
module load intel
module load impi
module load waitio

export WAITIO_MASTER_HOST=`hostname`-ib0
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=1

#waitio-serv-intel -m $WAITIO_MASTER_HOST
#echo  "serv host set ${WAITIO_MASTER_HOST} "

mpiexec -n 2  ./test-intel

exit
