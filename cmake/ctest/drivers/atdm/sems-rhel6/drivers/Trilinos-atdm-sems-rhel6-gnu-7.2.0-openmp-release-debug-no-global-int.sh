#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Experimental
fi
export ATDM_CONFIG_NO_GLOBAL_INT=ON
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/sems-rhel6/local-driver.sh
