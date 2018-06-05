#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Specialized
fi
export SALLOC_CTEST_TIME_LIMIT_MINUTES=0:15:00
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/chama/local-driver.sh
