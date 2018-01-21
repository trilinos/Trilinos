#!/bin/bash -l

/usr/bin/srun -N 1 --constraint=k80 -J $JOB_NAME \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/srun-job.sh

#   -e output.error
