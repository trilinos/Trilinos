#!/bin/bash
#export SLURM_CTEST_TIMEOUT=1:00:00
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=ATDM
fi
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/cts1/local-driver.sh
