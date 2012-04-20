#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on brain: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/usr/local/bin/git
export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Nightly

# CTest Settings
#

# Trilinos settings:
#

# Machine specific environment:
#

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on brain: `date`"
echo
