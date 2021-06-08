#!/bin/bash -l

if [ "${SLURM_CTEST_TIMEOUT}" == "" ] ; then
  SLURM_CTEST_TIMEOUT=1:20:00
  # This is just running tests, not the entire build!
fi

set -x

source $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-config-build.sh

set -x

atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver-test.sh \
  $PWD/ctest-s-driver-test.out \
  ${SLURM_CTEST_TIMEOUT}
