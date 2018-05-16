#!/bin/bash
export SRUN_CTEST_TIME_LIMIT_MINUTES=20
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/toss3/local-driver.sh
