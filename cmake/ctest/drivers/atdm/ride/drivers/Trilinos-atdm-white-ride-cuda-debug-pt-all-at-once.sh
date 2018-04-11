#!/bin/bash
#export Trilinos_TRACK=ATDM
export ATDM_CONFIG_ENABLE_ALL_PACKAGES=TRUE
export Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF
export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
export Trilinos_CTEST_USE_NEW_AAO_FEATURES=ON
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
