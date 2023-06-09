#!/bin/bash

#export SALLOC_CTEST_TIME_LIMIT_MINUTES=1:00:00

if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=PrimaryATDM
fi

# Disabling tests, see #10355
export CTEST_DO_TEST=OFF

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/van1-tx2/local-driver.sh
