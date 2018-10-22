#!/bin/bash

#
# This script only runs tests and it requires that the env to already to be
# set!  See README.md for detals.
#

set +x

echo
echo "Start: ctest-s-driver-test.sh"
echo
echo "  ==> `date`"
echo

echo "Running ctest -S comamnd to test only ..."

CTEST_S_CMND=env CTEST_DO_NEW_START=OFF CTEST_DO_UPDATES=OFF CTEST_DO_CONFIGURE=OFF CTEST_DO_BUILD=OFF CTEST_DO_TEST=ON ctest -V -S $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake

echo
echo "Running: ${CTEST_S_CMND} ..."

$CTEST_S_CMND
ATDM_TCD_CTEST_S_RETURN_CODE=$?

echo
echo "The 'ctest -S ctest-drivers.cmake' command returned code '$ATDM_TCD_CTEST_S_RETURN_CODE'"
echo
echo "End: ctest-s-driver-test.sh"
echo
echo "  ==> `date`"
echo
echo

# NOTE: The above output is important in order to determine if this script as
# run by the batch scheduler for the machine lets this script complete without
# killing it!
