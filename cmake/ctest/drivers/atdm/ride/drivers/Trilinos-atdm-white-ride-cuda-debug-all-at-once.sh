#!/bin/bash
#export Trilinos_TRACK=ATDM
export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
export Trilinos_CTEST_USE_NEW_AAO_FEATURES=ON
$WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ride/local-driver.sh
