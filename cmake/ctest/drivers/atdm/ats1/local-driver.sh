#!/bin/bash -l

if [ "${SBATCH_BUILD_TIME_LIMIT_MINUTES}" == "" ] ; then
  export SBATCH_BUILD_TIME_LIMIT_MINUTES=600 # Default 10 hour time limit
fi

if [ "${SBATCH_TEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  export SBATCH_TEST_TIME_LIMIT_MINUTES=840 # Default 14 hour time limit
fi

# Load environment on the login node
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Run configure, build on a haswell compute node
atdm_config_sbatch_extra_args="$ATDM_CONFIG_SBATCH_EXTRA_ARGS"
unset ATDM_CONFIG_SBATCH_EXTRA_ARGS
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver-build-on-allocation.sh \
  $PWD/ctest-s-driver.out \
  ${SBATCH_BUILD_TIME_LIMIT_MINUTES}

# Run tests on either a KNL or HSW compute node
# Per feedback in ATDV-347, disable running tests for now.
#export ATDM_CONFIG_SBATCH_EXTRA_ARGS="$atdm_config_sbatch_extra_args"
#echo "ATDM_CONFIG_SBATCH_EXTRA_ARGS = '${ATDM_CONFIG_SBATCH_EXTRA_ARGS}'"
#atdm_run_script_on_compute_node \
#  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver-test-on-allocation.sh \
#  $PWD/ctest-s-driver.out \
#  ${SBATCH_TEST_TIME_LIMIT_MINUTES} 
