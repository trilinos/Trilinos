#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Experimental
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/sems-rhel6/local-driver.sh
