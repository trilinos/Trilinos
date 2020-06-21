#!/bin/bash
#
# This script is designed to run inside of an allocation which then runs the
# build and tests on the same allocated compute node.  (See
# van1-tx2/local-drivers.sh and atdm/READM.md for details.)
#

set +x

echo
echo "======================================================================"
echo ""
echo "  Running ${JOB_NAME}"
echo "  in salloc allocation from node '$(hostname)'"
echo ""
echo "======================================================================"
echo

# Must use same site name for build and test results so they match up on
# CDash!
export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}

echo
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "***"
echo "*** Running start, update, configure and build on a compute node"
echo "***"
echo

set -x
env \
  CTEST_DO_TEST=FALSE \
  srun -N 1 \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
set +x

echo
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "***"
echo "*** Running tests on compute node launched from the login node"
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
