#!/bin/bash
export Trilinos_TRACK=ATDM
export SRUN_CTEST_TIME_LIMIT_MINUTES=480 # 8 hour time limit
export ATDM_CONFIG_USE_MAKEFILES=ON
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/local-driver.sh
