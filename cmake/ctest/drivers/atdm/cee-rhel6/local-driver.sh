#!/bin/bash -l

set +x

if [[ "${Trilinos_ENABLE_BUILD_STATS}" == "" ]] && \
   [[ ! $JOB_NAME == *"intel"* ]] \
  ; then
  export Trilinos_ENABLE_BUILD_STATS=ON
fi
echo "Trilinos_ENABLE_BUILD_STATS='${Trilinos_ENABLE_BUILD_STATS}'"

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
