#!/bin/bash
export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}
srun -N1  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
