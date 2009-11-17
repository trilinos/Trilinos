#!/bin/bash

# Source the 
#this is not the best thing to do, but we'll leave it for now
cd $HOME
source .bash_profile

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/Users/bmpersc/bin/ctest
EG_EXE=/Users/bmpersc/bin/eg
BASEDIR=/Users/bmpersc/nightly/Trilinos.base/development
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/s909348
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

echo
echo "Starting nightly Trilinos development testing on s909348: `date`"
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

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_opt_s909348.cmake -VV &> "MPI_OPT_DEV_$DATE.log"

echo
echo "Doing serial debug development build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_debug_s909348.cmake -VV &> "SERIAL_DEBUG_DEV_$DATE.log"

echo
echo "Ending nightly Trilinos development testing on s909348: `date`"
echo

