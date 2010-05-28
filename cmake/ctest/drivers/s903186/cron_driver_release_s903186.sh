#!/bin/bash

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/usr/local/CMake2.8.1/Contents/bin/ctest
EG_EXE=/Users/jmwille/bin/eg
BASEDIR=/Users/jmwille/TrilinosTestHarness/Trilinos10.2
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/s903186
BRANCH="trilinos-release-10-2-branch"
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"
export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

echo
echo "Starting nightly Trilinos Release 10.2 testing on s903186: `date`"
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
echo "Doing mpi optimized Release 10.2 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_opt_release_s903186.cmake -VV &> "MPI_OPT_RELEASE_10.2_$DATE.log"

echo
echo "Doing serial debug Release 10.2 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_debug_release_s903186.cmake -VV &> "SERIAL_DEBUG_RELEASE_10.2_$DATE.log"

echo
echo "Ending nightly Trilinos Release 10.2 testing on s903186: `date`"
echo

