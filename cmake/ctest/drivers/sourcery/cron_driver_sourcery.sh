#!/bin/bash

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/usr/bin/ctest
GIT_EXE=/usr/local/git/bin/git
BASEDIR=/Users/knmorri/NightlyTestingTrilinos/Trilinos
DRIVER_SCRIPT_DIR=$BASEDIR/cmake/ctest/drivers/sourcery
BRANCH="master"
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"
export PYTHONPATH=/usr/bin/python

echo
echo "Starting nightly Trilinos development testing on sourcery: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo

cd $BASEDIR
if [ -d Trilinos ]; then
  echo Doing an update of existing directory
  cd Trilinos
  $EG_EXE pull
  cd ..
else
  echo Cloning the repository because none exists yets
  $EG_EXE clone $TRILINOS_REPOSITORY_LOCATION
fi


echo
echo "Doing mpi optimized development build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_opt_nag_sourcery.cmake -VV &> "MPI_OPT_DEV_NAG_$DATE.log"

echo
echo "Doing serial debug development build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_debug_nag_sourcery.cmake -VV &> "SERIAL_DEBUG_DEV_NAG_$DATE.log"

echo
echo "Ending nightly Trilinos development testing on sourcery: `date`"
echo

