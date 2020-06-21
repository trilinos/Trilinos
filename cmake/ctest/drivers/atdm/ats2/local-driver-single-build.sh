#!/bin/bash
#
# This script is designed to run on the login node.  The start, update,
# configure and build are done on the login node.  This is because it is just
# not robust to try to run this on the compute node (see ATDV-324).  However,
# ctest running the tests themselves must be run from the launch node
# (allocated using bsub).  This is because you can only invoke jsrun from the
# launch node, not the login or compute node!
#

set +x

if [ "${LSF_CTEST_TEST_TIMEOUT}" == "" ] ; then
  export LSF_CTEST_TEST_TIMEOUT=4:00
fi

# Must use same site name for build and test results so they match up on
# CDash!
export CTEST_SITE=${ATDM_CONFIG_CDASH_HOSTNAME}

echo
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "***"
echo "*** Running [start], [update], configure and build on the login node"
echo "***"
echo

set -x
env \
  CTEST_DO_TEST=FALSE \
  $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
set +x

echo
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "***"
echo "*** Running tests from the launch node"
echo "***"
echo

set -x
env \
  CTEST_DO_NEW_START=FALSE \
  CTEST_DO_UPDATES=FALSE \
  CTEST_DO_CONFIGURE=FALSE \
  CTEST_DO_BUILD=FALSE \
  bsub -J ${ATDM_CONFIG_BUILD_NAME} -W ${LSF_CTEST_TEST_TIMEOUT} -Is \
    $WORKSPACE/Trilinos/cmake/ctest/drivers/atdm/ctest-s-driver.sh
set +x
