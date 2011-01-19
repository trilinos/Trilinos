#!/bin/bash

# This script must be run from this directory!

echo
echo "Starting nightly Trilinos testing on billmp1: `date`"
echo

SCRIPT_DIR=$PWD
cd $HOME
source $SCRIPT_DIR/bash_profile
cd -

#export TDD_PARALLEL_LEVEL=1
#export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"
#export TDD_CTEST_TEST_TYPE=Nightly

export TDD_PARALLEL_LEVEL=1
export TDD_CTEST_TEST_TYPE=Experimental

time /usr/bin/python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on billmp1: `date`"
echo

echo "Finished nightly Trilinos CMake tests billmp1: http://casl-dev.ornl.gov/cdash/index.php?project=Trilinos" | mailx -s "Nightly CTest: billmp1" bakercg@ornl.gov
