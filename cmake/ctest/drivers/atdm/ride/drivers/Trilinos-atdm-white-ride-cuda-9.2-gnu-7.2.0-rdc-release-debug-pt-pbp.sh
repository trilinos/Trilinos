#!/bin/bash
if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Specialized
fi
export ATDM_CONFIG_ENABLE_ALL_PACKAGES=TRUE
#export Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF
export ATDM_CONFIG_CONFIGURE_OPTIONS_FILES=cmake/std/atdm/ATDMDevEnvAllPtPackages.cmake
export Trilinos_CTEST_DO_ALL_AT_ONCE=FALSE
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
