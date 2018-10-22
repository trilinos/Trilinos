#!/bin/bash

#
# Direct Jenkins driver for GCC 4.8.3 + OpenMPI 1.10.1 + OpenMP build
#

DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

export CTEST_BUILD_CONFIGURATION_NAME=GCC-4.8.4-OpenMPI-1.10.1-MpiReleaseDebugSharedPtOpenMP

if [ "${CTEST_TEST_TYPE}" == "" ] ; then
  export CTEST_TEST_TYPE=Nighlty
fi

if [ "${Trilinos_TRACK}" == "" ] ; then
  export Trilinos_TRACK=Specialized
fi

${DRIVER_SCRIPT_DIR}/ctest_std_driver.sh
