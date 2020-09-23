#!/bin/bash

if [[ "${SALLOC_CTEST_TIME_LIMIT_MINUTES}" == "" ]] ; then
  export SALLOC_CTEST_TIME_LIMIT_MINUTES=3:00:00
fi

if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/tlcc2/local-driver.sh
