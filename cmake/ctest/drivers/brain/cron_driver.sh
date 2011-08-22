#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on brain: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/usr/local/git/bin/git
export TDD_PARALLEL_LEVEL=8
export TDD_CTEST_TEST_TYPE=Nightly
export TDD_IN_TESTING_MODE=ON

# CTest Settings
#
export CTEST_DO_SUBMIT=FALSE
export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=TRUE

# Trilinos settings:
#
export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on brain: `date`"
echo
