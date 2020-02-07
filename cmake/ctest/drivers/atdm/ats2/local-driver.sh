#!/bin/bash -l

if [ "${LSF_CTEST_TIMEOUT}" == "" ] ; then
  LSF_CTEST_TIMEOUT=4:00
  # This is just running tests, not the entire build!
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh

set -x

atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh \
  $PWD/ctest-s-driver-test.out \
  ${LSF_CTEST_TIMEOUT}

if [ "${Trilinos_CTEST_RUN_CUDA_AWARE_MPI}" == "1" ]; then
  export TPETRA_ASSUME_CUDA_AWARE_MPI=1
  export CTEST_BUILD_NAME=${ATDM_CONFIG_BUILD_NAME}_cuda-aware-mpi
  export CTEST_DO_UPDATES=FALSE
  #export CTEST_DO_NEW_START=FALSE
  export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
  export CTEST_DO_CONFIGURE=FALSE
  export CTEST_DO_BUILD=FALSE
  atdm_run_script_on_compute_node \
    $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh \
    $PWD/ctest-s-driver-test.out
    ${LSF_CTEST_TIMEOUT}
fi
