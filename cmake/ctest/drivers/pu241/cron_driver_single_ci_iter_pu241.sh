#!/bin/bash

echo
echo "Starting continuous integration Trilinos testing iteration on pu241: `date`"
echo

BASEDIR=/home/casl-vri-admin/CIDashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

# mpi compiler must have ifort in path
. $DRIVER_SCRIPT_DIR/load_intel_12.0.4_env

umask u=rwx,g=rwx,o=

export TDD_PARALLEL_LEVEL=4
export TDD_CTEST_TEST_TYPE=Continuous

# Submit the outer TDD tests to casl-dev always since these are CASL machines
export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/CDash/submit.php?project=TrilinosDriver"

export RUN_CI_SERVER=ON

time env python ../cron_driver.py

echo
echo "Ending continuous integration Trilinos testing iteration on pu241: `date`"
echo
