#!/bin/bash -l

set +x

if [[ "${Trilinos_REPOSITORY_LOCATION}" == "" ]] ; then
  export Trilinos_REPOSITORY_LOCATION=git@github.com:trilinos/Trilinos.git
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
