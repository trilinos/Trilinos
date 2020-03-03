#!/bin/bash

if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi

if [ "${SALLOC_CTEST_TIME_LIMIT_MINUTES}" == "" ] ; then
  export SALLOC_CTEST_TIME_LIMIT_MINUTES=540 # 9 hour time limit
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/mutrino/local-driver.sh
