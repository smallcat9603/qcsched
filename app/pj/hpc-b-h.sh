#!/bin/sh

#------ pjsub option --------#
#PJM -L rscgrp=qc-hpc-b-p-o
#PJM -L node=2
#PJM --mpi proc=4
#PJM -L elapse=1:00:00
#PJM -g gz00
#PJM -j

#------- Program execution -------#
mpiexec /work/gz00/z30130/qcsched/app/bin/lu.o 1024
