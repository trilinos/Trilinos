#!/bin/bash -l

if [ "${SALLOC_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SALLOC_CTEST_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh

set -x

/usr/local/bin/salloc -N 1 -p standard -J $JOB_NAME \
  --time=${SALLOC_CTEST_TIME_LIMIT_MINUTES} \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh
