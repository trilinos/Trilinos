#!/bin/csh

echo
echo "Starting nightly Trilinos development testing on dorksaber: `date`"
echo

#
# TrilinosDriver settings:
#

setenv TDD_PARALLEL_LEVEL 2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
setenv TDD_CTEST_TEST_TYPE Nightly

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#

setenv TDD_HTTP_PROXY "http://wwwproxy.sandia.gov:80"
setenv http_proxy "http://wwwproxy.sandia.gov:80"
setenv TDD_FORCE_CMAKE_INSTALL 1

source ~/.cshrc

module load gcc/4.7.2

# Get a new git further up in the path
setenv PATH "/usr/local/bin":"$PATH"

# Matlab goodies
setenv LD_PRELOAD "/usr/local/gcc/gcc-4.7.2/lib64/libstdc++.so.6 /usr/local/gcc/gcc-4.7.2/lib64/libgcc_s.so.1 /usr/local/gcc/gcc-4.7.2/lib64/libgfortran.so.3"

env

which git

# Machine independent cron_driver:
#

setenv SCRIPT_DIR `dirname "$0"`
echo "SCRIPT_DIR = " $SCRIPT_DIR
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on dorksaber: `date`"
echo
