#!/bin/bash -l

set +x

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi

if [[ "${Trilinos_REPOSITORY_LOCATION}" == "" ]] ; then
  export Trilinos_REPOSITORY_LOCATION=git@github.com:trilinos/Trilinos.git
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
