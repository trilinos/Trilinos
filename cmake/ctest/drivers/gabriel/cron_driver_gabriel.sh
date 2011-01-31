#!/bin/bash

# This script must be run from this directory!

echo
echo "Starting nightly Trilinos testing on gabriel: `date`"
echo

SCRIPT_DIR=$PWD
cd $HOME
source $SCRIPT_DIR/bash_profile
cd -

export TDD_PARALLEL_LEVEL=2
#export TDD_PARALLEL_LEVEL=1
export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80/"
export TDD_CTEST_TEST_TYPE=Nightly
time /usr/local/bin/python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on gabriel: `date`"
echo

/home/rabartl/mailmsg.py "Finished nightly Trilinos CMake tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
