#!/bin/sh

#------ pjsub option --------#
#PJM -L rscgrp=qc-hpc-ibm-a
#PJM -L gpu=1
#PJM --mpi proc=4
#PJM -L elapse=1:00:00
#PJM -g gz00
#PJM -j

#------- Program execution -------#
module load gcc ompi
mpiexec -machinefile $PJM_O_NODEINF -n $PJM_MPI_PROC -npernode 4 /work/gz00/z30130/qcsched/app/bin/lu.a 1024
