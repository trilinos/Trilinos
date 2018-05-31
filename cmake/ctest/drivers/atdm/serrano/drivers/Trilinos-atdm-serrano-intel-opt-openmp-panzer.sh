#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
export SALLOC_CTEST_TIME_LIMIT_MINUTES=0:30:00
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/serrano/local-driver.sh
