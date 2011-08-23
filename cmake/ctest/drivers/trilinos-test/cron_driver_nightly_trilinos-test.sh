#!/bin/bash

CTEST_EXE=/home/trilinos/cmake/bin/ctest
EG_EXE=/home/trilinos/git/bin/eg
BASEDIR=/home/bmpersc/nightly/Trilinos.base/release_10
BASEDATADIR=/home/bmpersc/nightly
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test
BRANCH="trilinos-release-10-0-branch"
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

export TRILINOSDATADIRECTORY=$BASEDATADIR/TrilinosData

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

echo
echo "Starting nightly Trilinos Release 10.0 testing on trilinos-test: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo


cd $BASEDATADIR
#checkout the trilinos data directory. Needed for a few tests to run
cvs -q -d :ext:software:/space/CVS co TrilinosData

cd $BASEDIR
#checkout the bits of trilinos needed for running the nightly test scripts
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

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_release_trilinos-test.cmake -VV &> "MPI_RELEASE_10.0_$DATE.log"

echo
echo "Doing mpi optimized release 10.0  shared library build: `date`"
echo 

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_release_shared_trilinos-test.cmake -VV &> "MPI_RELEASE_10.0_SHARED_$DATE.log"

echo
echo "Doing serial release 10.0 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_release_trilinos-test.cmake -VV &> "SERIAL_RELEASE_10.0_$DATE.log"

echo
echo "Ending nightly Trilinos testing on trilinos-test: `date`"
echo
