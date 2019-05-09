#!/bin/bash
export SALLOC_CTEST_TIME_LIMIT_MINUTES=0:45:00
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/serrano/local-driver.sh
