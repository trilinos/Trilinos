#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Experimental
fi
export SBATCH_TIME_LIMIT_MINUTES=900 # 15 hour time limit
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver.sh
