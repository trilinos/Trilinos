#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=PrimaryATDM
fi

$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats2/local-driver.sh
