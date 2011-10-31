#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on ictinus: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/Users/jmwille/install/git/bin/git
export TDD_PARALLEL_LEVEL=1

# Trilinos settings:
#

#export CTEST_TEST_TYPE=Experimental

#export CTEST_DO_SUBMIT=FALSE

#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

export PATH=/Users/jmwille/install/git/bin:/Users/jmwille/install/cmake-2.8.4/bin:/Users/jmwille/install/gcc-4.4.6/openmpi-1.4.3/bin:/Users/jmwille/install/gcc-4.4.6/bin:$PATH

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on ictinus: `date`"
echo
