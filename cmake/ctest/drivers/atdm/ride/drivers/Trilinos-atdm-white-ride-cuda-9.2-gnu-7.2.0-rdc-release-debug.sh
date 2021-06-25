#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=PrimaryATDM
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
