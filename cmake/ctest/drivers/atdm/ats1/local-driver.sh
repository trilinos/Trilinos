#!/bin/bash -l

if [ "${SBATCH_TIME_LIMIT_MINUTES}" == "" ] ; then
  SBATCH_TIME_LIMIT_MINUTES=180  # Default limit is 3 hours
fi

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

# Load environment on the login node
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

set -x

#export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}

# Run configure, build, and test on the compute node
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver-on-allocation.sh \
  $PWD/ctest-s-driver.out \
  ${SBATCH_TIME_LIMIT_MINUTES}
