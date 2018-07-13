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


# Machine specific environment
#

setenv TDD_HTTP_PROXY "http://wwwproxy.sandia.gov:80"
setenv http_proxy "http://wwwproxy.sandia.gov:80"
setenv TDD_FORCE_CMAKE_INSTALL 1
setenv TDD_DEBUG_VERBOSE 1

source ~/.cshrc

# Get a new git further up in the path
setenv PATH "/usr/local/bin":"$PATH"

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
