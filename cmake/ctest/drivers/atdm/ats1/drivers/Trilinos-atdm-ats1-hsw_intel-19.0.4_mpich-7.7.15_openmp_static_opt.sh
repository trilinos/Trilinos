#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=SparcATDM
fi
export Trilinos_SKIP_CTEST_ADD_TEST=FALSE
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ats1/local-driver.sh
