#!/bin/bash
export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git
export SRUN_CTEST_TIME_LIMIT_MINUTES=300 # 5 hour time limit
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/shiller/jenkins-driver.sh
