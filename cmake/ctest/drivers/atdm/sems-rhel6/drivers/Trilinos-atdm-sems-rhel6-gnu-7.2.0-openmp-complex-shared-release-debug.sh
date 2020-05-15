#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/sems-rhel6/local-driver.sh
