#!/bin/bash -l

if [ "${LSF_CTEST_ALL_TIMEOUT}" == "" ] ; then
  LSF_CTEST_ALL_TIMEOUT=12:00
fi

if [ "${LSF_CTEST_TEST_TIMEOUT}" == "" ] ; then
  LSF_CTEST_TEST_TIMEOUT=4:00
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

set -x

atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh \
  $PWD/ctest-s-driver.out \
  ${LSF_CTEST_ALL_TIMEOUT}

if [ "${Trilinos_CTEST_RUN_CUDA_AWARE_MPI}" == "1" ]; then
  # Wait a little while before attempting another bsub allocation
  # to try working around "'Error: Remote JSM server is not responding'"
  sleep 5
  export CTEST_BUILD_NAME=${ATDM_CONFIG_BUILD_NAME}_cuda-aware-mpi
  export TPETRA_ASSUME_CUDA_AWARE_MPI=1
  export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
  export CTEST_DO_UPDATES=OFF
  export CTEST_DO_CONFIGURE=OFF
  export CTEST_DO_BUILD=OFF
  atdm_run_script_on_compute_node \
    $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh \
    $PWD/ctest-s-driver-cuda-aware-mpi.out \
    ${LSF_CTEST_TEST_TIMEOUT}
fi
