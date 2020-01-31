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
