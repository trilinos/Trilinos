#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on typhon: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#

export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"
export TDD_FORCE_CMAKE_INSTALL=1
export CUDA_LAUNCH_BLOCKING=1
export OMP_NUM_THREADS=2
export PATH="$PATH:/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/intel64"
# Machine independent cron_driver:
#
#openmpi-1.7-cuda6

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on typhon: `date`"
echo
