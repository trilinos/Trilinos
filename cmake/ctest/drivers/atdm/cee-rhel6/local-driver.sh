#!/bin/bash -l

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
