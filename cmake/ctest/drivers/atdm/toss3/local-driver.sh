#!/bin/bash -l

if [ "${SRUN_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SRUN_CTEST_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

set -x

export CTEST_DO_TEST=OFF
#/usr/bin/srun -N 1 --time=600 --account=fy150090 -J $JOB_NAME \
#  --time=${SRUN_CTEST_TIME_LIMIT_MINUTES} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
