#!/bin/bash
if [[ "${CTEST_DO_TEST}" == "" ]] ; then
  export CTEST_DO_TEST=OFF
fi

export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}
srun -N1  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
