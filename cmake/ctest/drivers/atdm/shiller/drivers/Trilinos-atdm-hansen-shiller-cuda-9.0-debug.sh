#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
export SRUN_CTEST_TIME_LIMIT_MINUTES=480 # 8 hour time limit
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/local-driver.sh
