#!/bin/bash -l

if [ "${SRUN_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SRUN_CTEST_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

set -x

export CTEST_DO_TEST=OFF

#/usr/local/bin/salloc -N 1 -p standard -J $JOB_NAME \
#  --time=${SRUN_CTEST_TIME_LIMIT_MINUTES} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
