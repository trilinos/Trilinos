#!/bin/bash

echo
echo "Starting nightly Trilinos testing on pu241: `date`"
echo

BASEDIR=/home/casl-vri-admin/TrilinosNightlyDashboardsIntel12
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/pu241
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"


# Allow override of dev env
if [ "$CASL_VRI_DEV_ENV_BASE" == "" ] ; then
  CASL_VRI_DEV_ENV_BASE=/opt/casl_vri_dev_env
fi
echo "CASL_VRI_DEV_ENV_BASE = '$CASL_VRI_DEV_ENV_BASE'"

# load intel 12 dev env
. $CASL_VRI_DEV_ENV_BASE/fissile_four/build_scripts/load_official_dev_env.sh 12
# mark environment so that dashboard will build the intel 12 configs only
export INTEL12_BUILD=1

umask u=rwx,g=rwx,o=

export TDD_PARALLEL_LEVEL=4
export TDD_CTEST_TEST_TYPE=Nightly

# Submit the outer TDD tests to casl-dev always since these are CASL machines
export TDD_CTEST_DROP_SITE=casl-dev.ornl.gov
export TDD_CTEST_DROP_LOCATION="/cdash/submit.php?project=TrilinosDriver"

export TDD_BUILD_NAME="Linux-TDD-Intel12-pu241"
#export CTEST_TEST_TYPE=Experimental

time env python ../cron_driver.py

echo
echo "Ending nightly Trilinos testing on pu241 (Intel 12): `date`"
echo

