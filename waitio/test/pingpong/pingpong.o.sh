#!/bin/sh

#------ pjsub option --------#
#PJM -L rscgrp=regular-o
#PJM -L node=2
#PJM --mpi proc=2
#PJM -L elapse=1:00:00
#PJM -g gz00
#PJM -j

#------- Program execution -------#
mpiexec /work/gz00/z30130/mpi/waitio/test/pingpong/pingpong
