#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=PrimaryATDM
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/sems-rhel7/local-driver.sh
