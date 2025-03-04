#!/bin/sh

#------ pjsub option --------#
#PJM -L rscgrp=regular-a
#PJM -L node=2
#PJM --mpi proc=2
#PJM -L elapse=1:00:00
#PJM -g gz00
#PJM -j

#------- Program execution -------#
module load gcc ompi
mpiexec -machinefile $PJM_O_NODEINF -n $PJM_MPI_PROC -npernode 1 /work/gz00/z30130/mpi/waitio/test/pingpong/pingpong
