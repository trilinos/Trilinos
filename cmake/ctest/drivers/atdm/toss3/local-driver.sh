#!/bin/bash -l

if [ "${SALLOC_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SALLOC_CTEST_TIME_LIMIT_MINUTES=1:30:00
  # This is just running tests, not the entire build!
fi

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh

set -x

salloc -N1 --time=${SALLOC_CTEST_TIME_LIMIT_MINUTES} --account=fy150090 -J $JOB_NAME \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh
