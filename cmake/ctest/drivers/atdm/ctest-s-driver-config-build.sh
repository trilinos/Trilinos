#!/bin/bash -l

set +x

echo
echo "Start: ctest-s-driver-config-build.sh"
echo
echo "  ==> `date`"
echo

echo "Loading env and running ctest -S comamnd to update, configure and build ..."

source ${WORKSPACE}/Trilinos/cmake/ctest/drivers/atdm/utils/setup_env.sh

source ${WORKSPACE}/Trilinos/cmake/ctest/drivers/atdm/utils/create-src-and-build-dir.sh

CTEST_S_CMND="env CTEST_DO_TEST=OFF ctest -V -S $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-driver.cmake"

echo
echo "Running: \"${CTEST_S_CMND}\" ..."

$CTEST_S_CMND
ATDM_TCD_CTEST_S_RETURN_CODE=$?

echo
echo "The 'ctest -S ctest-drivers.cmake' command returned code '$ATDM_TCD_CTEST_S_RETURN_CODE'"
echo
echo "End: ctest-s-driver-config-build.sh"
echo
echo "  ==> `date`"
echo
echo

# NOTE: The above output is important in order to determine if this script as
# run by the batch scheduler for the machine lets this script complete without
# killing it!
