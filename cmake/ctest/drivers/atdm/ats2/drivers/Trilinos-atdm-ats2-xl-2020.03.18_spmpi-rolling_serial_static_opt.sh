#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver.sh
