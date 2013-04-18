#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on eleanor: `date`"
echo

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

export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80"
export http_proxy="http://wwwproxy.sandia.gov:80"
export TDD_FORCE_CMAKE_INSTALL=0

source /usr/local/intel/Compiler/11.1/064/bin/iccvars.sh intel64
source /usr/local/intel/Compiler/11.1/064/mkl/tools/environment/mklvarsem64t.sh

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on eleanor: `date`"
echo
