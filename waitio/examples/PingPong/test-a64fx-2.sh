#!/bin/sh
#------ pjsub option --------#
#PJM -N "test1"
#PJM -L rscgrp=coupler-o
#PJM -L node=1
#PJM --mpi proc=2
#PJM -L elapse=00:02:00
#PJM -g gz00
#PJM -j
#------- Program execution -------#

hostname
module purge
module load fj
module load fjmpi
module load waitio

export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=1
export WAITIO_NPB=2
export WAITIO_MASTER_HOST=`waitio-serv-a64fx -c`

echo  "serv host is ${WAITIO_MASTER_HOST} "
mpiexec  ./test-a64fx
exit
