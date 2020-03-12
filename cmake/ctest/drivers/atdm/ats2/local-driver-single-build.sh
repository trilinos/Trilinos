#!/bin/bash
#
# This script is designed to run on the launch node using a single node
# allocation.
#

set +x

# Must get the same site name for both builds show they show up on same CDash build
export CTEST_SITE=$(echo ${LSB_HOSTS} | cut -d' ' -f2)

echo
echo "***"
echo "*** Running [start], [update], configure and build for job '${LSB_JOBNAME}'"
echo "***"
echo

set -x
env \
  CTEST_DO_TEST=FALSE \
  lrun -n 1 \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
set +x

echo
echo "***"
echo "*** Running tests for job '${LSB_JOBNAME}'"
echo "***"
echo

set -x
env \
  CTEST_DO_NEW_START=FALSE \
  CTEST_DO_UPDATES=FALSE \
  CTEST_DO_CONFIGURE=FALSE \
  CTEST_DO_BUILD=FALSE \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
set +x
