#!/bin/bash

echo
echo "Starting nightly Trilinos testing on u233: `date`"
echo

source /opt/casldev/env/casl_dev_env.sh

BASEDIR=/home/casl-vri-admin/Dashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/u233
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Experimental

export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/CDash/submit.php?project=TrilinosDriver"

time env python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on u233: `date`"
echo

echo "Finished nightly Trilinos CMake tests u233: http://casl-dev.ornl.gov/CDash/index.php?project=Trilinos" | mailx -s "Nightly CTest: u233" bakercg@ornl.gov

# ToDo: Above: add rabartl@sandia.gov and others later?
