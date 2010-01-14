#!/bin/bash

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/usr/local/bin/ctest
EG_EXE=/Users/jmwille/bin/eg
BASEDIR=/Users/jmwille/TrilinosTestHarness
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/s903186
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

echo
echo "Starting nightly Trilinos development testing on s903186: `date`"
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

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_opt_s903186.cmake -VV &> "MPI_OPT_DEV_$DATE.log"

echo
echo "Doing serial debug development build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_debug_s903186.cmake -VV &> "SERIAL_DEBUG_DEV_$DATE.log"

echo
echo "Ending nightly Trilinos development testing on s903186: `date`"
echo

