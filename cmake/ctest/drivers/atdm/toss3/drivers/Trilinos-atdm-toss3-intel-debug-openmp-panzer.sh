#!/bin/bash
export SALLOC_CTEST_TIME_LIMIT_MINUTES=0:45:00
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/toss3/local-driver.sh
