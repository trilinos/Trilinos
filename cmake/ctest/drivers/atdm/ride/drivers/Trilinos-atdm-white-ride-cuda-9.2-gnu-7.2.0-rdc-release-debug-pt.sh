#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Specialized
fi
if [ "${ATDM_CONFIG_KNOWN_HOSTNAME}" == "ride" ] ; then
  export EXCLUDE_NODES_FROM_BSUB="-R hname!=ride7"
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
