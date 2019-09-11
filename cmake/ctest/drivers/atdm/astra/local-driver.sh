#!/bin/bash -l

if [ "${SALLOC_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  SALLOC_CTEST_TIME_LIMIT_MINUTES=1:30:00
  # This is just running tests, not the entire build!
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh

set -x



