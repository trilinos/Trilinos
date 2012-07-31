#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on sourcery: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/usr/local/git/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_HTTP_PROXY="http://wwwproxy.ca.sandia.gov:80/"

# Trilinos settings:
#

export TDD_CTEST_TEST_TYPE=Nightly
export CTEST_TEST_TYPE=Experimental
export CTEST_DO_SUBMIT=ON
export TDD_DO_SUBMIT=ON
export Trilinos_PACKAGES=ForTrilinos 
export CTEST_DO_UPDATES=TRUE
export TDD_IN_TESTING_MODE=ON 
export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE
#export TDD_FORCE_CMAKE_INSTALL=0 
#export TDD_FORCE_INNER_CMAKE_INSTALL=0 
export CTEST_EXPLICITLY_ENABLE_IMPLICITLY_ENABLED_PACKAGES=TRUE 

# Machine specific environment:
#

export PYTHONPATH=/usr/bin/python

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on sourcery: `date`"
echo
