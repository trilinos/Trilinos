#!/bin/bash -l

if [ "${SRUN_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SRUN_CTEST_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

/usr/bin/srun -N 1 --constraint=k80 -J $JOB_NAME \
  --time=${SRUN_CTEST_TIME_LIMIT_MINUTES} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
