#!/bin/bash
export SALLOC_CTEST_TIME_LIMIT_MINUTES=0:45:00
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Specialized
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/chama/local-driver.sh
