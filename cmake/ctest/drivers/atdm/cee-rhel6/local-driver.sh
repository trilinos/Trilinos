#!/bin/bash -l

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
