#!/bin/bash -l

if [ "${BSUB_CTEST_TIME_LIMIT}" == "" ] ; then
  export BSUB_CTEST_TIME_LIMIT=04:00:00
fi

set -x

bsub -x -I \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/bsub-ctest-job.sh
