#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Experimental
fi
export TPETRA_ASSUME_CUDA_AWARE_MPI=0
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver.sh

export TPETRA_ASSUME_CUDA_AWARE_MPI=1
export CTEST_BUILD_NAME=Trilinos-atdm-ats2-cuda-10.1.243-gnu-7.3.1-spmpi-2019.06.24_cuda-aware_static_dbg
export CTEST_DO_UPDATES=FALSE
#export CTEST_DO_NEW_START=FALSE
export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
export CTEST_DO_CONFIGURE=FALSE
export CTEST_DO_BUILD=FALSE
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh \
  $PWD/ctest-s-driver-test.out
