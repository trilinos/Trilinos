#!/bin/bash

echo
echo "Starting nightly Trilinos testing on pu241: `date`"
echo

BASEDIR=/home/casl-vri-admin/TrilinosNightlyDashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

. /projects/vera/gcc-4.6.1/load_dev_env.sh

umask u=rwx,g=rwx,o=

export TDD_PARALLEL_LEVEL=6
export TDD_CTEST_TEST_TYPE=Nightly

# Submit the outer TDD tests to casl-dev always since these are CASL machines
export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/cdash/submit.php?project=TrilinosDriver"

#export CTEST_TEST_TYPE=Experimental

time env python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on pu241: `date`"
echo

