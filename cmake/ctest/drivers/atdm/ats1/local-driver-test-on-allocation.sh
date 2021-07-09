#!/bin/bash

echo
echo "======================================================================="
echo "***"
echo "*** Run tests on compute node '${SLURM_JOB_NODELIST}'"
echo "***"
echo

export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}
export CTEST_DO_NEW_START=FALSE
export CTEST_DO_UPDATES=FALSE
export CTEST_DO_CONFIGURE=FALSE
export CTEST_DO_BUILD=FALSE
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
