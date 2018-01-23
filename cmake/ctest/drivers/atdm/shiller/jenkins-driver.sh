#!/bin/bash -l

if [ "${SRUN_TIME_LIMIT_MIN}" == "" ] ; then
  SRUN_TIME_LIMIT_MIN=180  # Default limit is 3 hours
fi

set -x

/usr/bin/srun -N 1 --constraint=k80 -J $JOB_NAME \
  --time-min=${SRUN_TIME_LIMIT_MIN} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/srun-job.sh

#   -e output.error
