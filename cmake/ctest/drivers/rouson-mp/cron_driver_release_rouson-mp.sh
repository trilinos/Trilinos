#!/bin/bash

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/usr/bin/ctest
GIT_EXE=/usr/local/git/bin/git
BASEDIR=/Users/knmorri/NightlyTestingTrilinos/Trilinos
DRIVER_SCRIPT_DIR=$BASEDIR/cmake/ctest/drivers/rouson-mp
BRANCH="trilinos-release-10-10-branch"
#BRANCH="master"
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"
export PYTHONPATH=/usr/bin/python

echo
echo "Starting nightly Trilinos Release 10.10 testing on rouson-mp: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo

cd $BASEDIR
if [ -d Trilinos ]; then
  echo Doing an update of existing directory
  cd Trilinos
  $GIT_EXE pull
  cd ..
else
  echo Cloning the repository because none exists yets
  $GIT_EXE clone $TRILINOS_REPOSITORY_LOCATION
fi


echo
echo "Doing mpi optimized Release 10.10 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_opt_release_nag_rouson-mp.cmake -VV &> "MPI_OPT_RELEASE_NAG_10.10_$DATE.log"

echo
echo "Doing serial debug Release 10.10 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_debug_release_nag_rouson-mp.cmake -VV &> "SERIAL_DEBUG_RELEASE_NAG_10.10_$DATE.log"

echo
echo "Ending nightly Trilinos Release 10.10 testing on rouson-mp: `date`"
echo

