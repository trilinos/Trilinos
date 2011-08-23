#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on godel: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/home/trilinos/install/bin/eg
export TDD_PARALLEL_LEVEL=4
if [ -z "$TDD_CTEST_TEST_TYPE" ]; then
    export TDD_CTEST_TEST_TYPE=Nightly
fi

# Trilinos settings:
#

if [ -z "$CTEST_TEST_TYPE" ]; then
    export CTEST_TEST_TYPE=Nightly
fi

#export CTEST_TEST_TYPE=Experimental

#export CTEST_DO_SUBMIT=FALSE

#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on godel: `date`"
echo
