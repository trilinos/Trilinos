#!/bin/bash
#
# This script is designed to run on the launch node using a single node
# allocation.  This script is called by the ats2/local-drivers.sh script in a
# bsub allocation as:
#
#  bsub -J <job-name> -W <timeout> -Is <this-dir>/local-driver-single-build.sh
#
# Therefore, this script runs on the launch node (not the login node or the
# compute node).
# 
# This script matches the approach documented in the cmake/std/atdm/README.md
# file in the "ATS-2" section where everything except the tests are run on the
# compute node using 'lrun', and then the tests themselves are launched
# directly from the launch node (since 'jsrun' can only be run from the launch
# node, not the login node or the compute node).
#

set +x

# Must get the same site name for both builds show they show up on same CDash build
if [[ "${CTEST_TEST_TYPE}" == "Experimental" ]] ; then
  export CTEST_SITE=$(atdm_ats2_get_allocated_compute_node_name)
  # NOTE: Must only conditionally set this because if set in the env, will
  # override the default value set in the file
  # TrilinosCTestDriverCore.atdm.cmake!
fi

echo
echo "***"
echo "*** Running [start], [update], configure and build for job '${LSB_JOBNAME}'"
echo "*** on the compute node"
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
echo "*** Running tests for job '${LSB_JOBNAME}' from the launch node"
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
