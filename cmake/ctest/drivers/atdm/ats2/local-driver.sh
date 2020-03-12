#!/bin/bash -l

if [ "${LSF_CTEST_ALL_TIMEOUT}" == "" ] ; then
  LSF_CTEST_ALL_TIMEOUT=12:00
fi

if [ "${LSF_CTEST_TEST_TIMEOUT}" == "" ] ; then
  LSF_CTEST_TEST_TIMEOUT=4:00
fi

# Need to load env so we define some vars
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Allow default setting for TPETRA_ASSUME_CUDA_AWARE_MPI=0 in trilinos_jsrun
unset TPETRA_ASSUME_CUDA_AWARE_MPI

set -x
bsub -J ${ATDM_CONFIG_BUILD_NAME} -W ${LSF_CTEST_ALL_TIMEOUT} -Is \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver-single-build.sh
set -x

if [[ "${Trilinos_CTEST_RUN_CUDA_AWARE_MPI}" == "1" ]]; then
  # Wait a little while before attempting another bsub allocation
  # to try working around "'Error: Remote JSM server is not responding'"
  sleep 5
  export CTEST_BUILD_NAME=${ATDM_CONFIG_BUILD_NAME}_cuda-aware-mpi
  if [[ "${CTEST_TEST_TYPE}" == "Experimental" ]] ; then
    export CTEST_BUILD_NAME="${CTEST_BUILD_NAME}-exp"
  fi
  export TPETRA_ASSUME_CUDA_AWARE_MPI=1
  export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE
  export CTEST_DO_UPDATES=OFF
  set -x
  bsub -J ${CTEST_BUILD_NAME} -W ${LSF_CTEST_TEST_TIMEOUT} -Is \
    $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver-single-build.sh
  set +x
fi

# NOTE: Above, we are using a single bsub allocation for each job that appears
# on CDash.  We run the configure, build, the tests within that allocation.
# We do this so that the hostname that shows up in the cmake configure output
# that is posted to CDash is the same node that runs the test.  However, we
# use a seprate node allocation to run the cuda-aware MPI build and tests.
# That is so that if we are hit by the jsrun bug in the first non-cuda-aware
# MPI test runs (see #6855), then this will not break the tests in the
# cuda-aware-mpi runs.  Otherwise, you could do a single bsub allocation to
# run both the non-cuda-aware MPI and the cuda-aware MPI builds.  But the
# jsrun bug is why we don't do that.
#
# NOTE: Above, we put an '-exp' on the end of the cuda-aware-mpi build name to
# make it more clear that this is an experimental build!
