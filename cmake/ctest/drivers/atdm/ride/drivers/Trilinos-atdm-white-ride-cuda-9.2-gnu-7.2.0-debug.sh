#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
export EXCLUDE_NODES_FROM_BSUB="-R hname!=ride14"
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
