#!/bin/sh

#------ pjsub option --------#
#PJM -L rscgrp=regular-o
#PJM -L node=2x4:torus
#PJM --mpi proc=16
#PJM -L elapse=1:00:00
#PJM -g gz00
#PJM -j

#------- Program execution -------#
mpiexec /work/gz00/z30130/qcsched/app/bin/lu.o 1024
