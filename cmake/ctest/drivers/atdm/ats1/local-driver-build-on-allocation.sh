#!/bin/bash

echo
echo "======================================================================="
echo "***"
echo "*** Run update, configure and build on compute node '${SLURM_JOB_NODELIST}'"
echo "***"
echo

export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}
export CTEST_DO_TEST=FALSE
srun -N1 $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
