#!/bin/bash -l

set +x

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
