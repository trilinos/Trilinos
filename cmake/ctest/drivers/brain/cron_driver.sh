#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on brain: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/usr/local/git/bin/git
export TDD_PARALLEL_LEVEL=8

# CTest Settings
#

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
