#!/bin/bash -l

set +x

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

if [[ "${Trilinos_SKIP_CTEST_ADD_TEST}" == "" ]] && \
   [[ "${ATDM_CONFIG_USE_CUDA}" == "ON" ]] \
  ; then
  # For CUDA builds, do not run the the test and example executables to speed
  # things up and to allow us to build and install on CEE machines that don't
  # actaully have a GPU!
  export Trilinos_SKIP_CTEST_ADD_TEST=ON
  # NOTE: If building the tests and examples is too expensive, we could set
  # Trilinos_INNER_ENABLE_TESTS=OFF.
fi

set -x

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
