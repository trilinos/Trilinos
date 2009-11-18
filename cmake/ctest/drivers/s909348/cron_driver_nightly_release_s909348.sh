#!/bin/bash

# Source the 
#this is not the best thing to do, but we'll leave it for now
cd $HOME
source .bash_profile

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/Users/bmpersc/bin/ctest
EG_EXE=/Users/bmpersc/bin/eg
BASEDIR=/Users/bmpersc/nightly/Trilinos.base/release_10
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/s909348
BRANCH="trilinos-release-10-0-branch"
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

echo
echo "Starting nightly Trilinos release 10.0 testing on s909348: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo


cd $BASEDIR
if [ -d Trilinos ]; then
  echo Doing an update of existing directory
  cd Trilinos
  $EG_EXE switch $BRANCH
  $EG_EXE pull
  cd ..
else
  echo Cloning the repository because none exists yets
  $EG_EXE clone $TRILINOS_REPOSITORY_LOCATION
  cd Trilinos
  $EG_EXE switch $BRANCH
  cd ..
fi

echo
echo "Doing mpi optimized release 10.0 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_mpi_release_s909348.cmake -VV &> "MPI_RELEASE_10.0_$DATE.log"

echo
echo "Doing serial release 10.0 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_mac_nightly_serial_release_s909348.cmake -VV &> "SERIAL_RELEASE_10.0_$DATE.log"

echo
echo "Ending nightly Trilinos release 10.0 testing on s909348: `date`"
echo

