#!/bin/bash

echo
echo "Starting nightly Trilinos testing on pu241: `date`"
echo

. /usr/share/Modules/init/bash
module use /opt/casl_vri_dev_env/fissile_four/modules
module load casl-dev
module load intel/12.0.4
umask u=rwx,g=rwx,o=

BASEDIR=/home/casl-vri-admin/Dashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

export TDD_PARALLEL_LEVEL=4
export TDD_CTEST_TEST_TYPE=Nightly

# Submit the outer TDD tests to casl-dev always since these are CASL machines
export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/CDash/submit.php?project=TrilinosDriver"

time env python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on pu241: `date`"
echo

echo "Finished nightly Trilinos CMake tests pu241: http://casl-dev.ornl.gov/CDash/index.php?project=Trilinos" | mailx -s "Nightly CTest: pu241" casl-vri-admin@localhost
