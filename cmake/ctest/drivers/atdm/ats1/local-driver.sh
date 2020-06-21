#!/bin/bash -l

if [ "${SBATCH_BUILD_TIME_LIMIT_MINUTES}" == "" ] ; then
  export SBATCH_BUILD_TIME_LIMIT_MINUTES=600 # Default 10 hour time limit
fi

if [ "${SBATCH_TEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  export SBATCH_TEST_TIME_LIMIT_MINUTES=780 # Default 13 hour time limit
fi

# Load environment on the login node
source $WORKSPACE/Trilinos/cmake/std/atdm/load-env.sh $JOB_NAME

# Don't run any tests for now until we can get them to run faster
if [[ "${Trilinos_SKIP_CTEST_ADD_TEST}" == "" ]] ; then
  export Trilinos_SKIP_CTEST_ADD_TEST=TRUE
fi

# Don't install the intel-18 builds because SPARC does not support intel-18 on
# this platform and EMPIRE does not use these installs (yet) (ATDV-361)
if [[ "${ATDM_CONFIG_COMPILER}" = "INTEL-18"* ]] ; then
  export CTEST_DO_INSTALL=OFF
fi

if [[ "${Trilinos_INNER_ENABLE_TESTS}" == "" ]] && \
   [[ "${ATDM_CONFIG_BUILD_TYPE}" == "DEBUG" ]] \
  ; then
  # For debug builds, not build the test and example executables to speed
  # things up for now.
  export Trilinos_INNER_ENABLE_TESTS=OFF
fi

# Run configure, build on a haswell compute node
atdm_config_sbatch_extra_args="$ATDM_CONFIG_SBATCH_EXTRA_ARGS"
unset ATDM_CONFIG_SBATCH_EXTRA_ARGS
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver-build-on-allocation.sh \
  $PWD/ctest-s-driver.out \
  ${SBATCH_BUILD_TIME_LIMIT_MINUTES}

# Run tests on either a KNL or HSW compute node
# Per feedback in ATDV-347, disable running tests for now.
export ATDM_CONFIG_SBATCH_EXTRA_ARGS="$atdm_config_sbatch_extra_args"
echo "ATDM_CONFIG_SBATCH_EXTRA_ARGS = '${ATDM_CONFIG_SBATCH_EXTRA_ARGS}'"
atdm_run_script_on_compute_node \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver-test-on-allocation.sh \
  $PWD/ctest-s-driver.out \
  ${SBATCH_TEST_TIME_LIMIT_MINUTES}
