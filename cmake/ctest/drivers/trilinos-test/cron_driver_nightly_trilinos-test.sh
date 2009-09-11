#!/bin/bash

CTEST_EXE=/home/trilinos/cmake/bin/ctest
BASEDIR=/home/bmpersc/nightly/Trilinos.base
CHECKOUT_DIR=Trilinos_release_10.0
DRIVER_SCRIPT_DIR=$BASEDIR/$CHECKOUT_DIR/cmake/ctest/drivers/trilinos-test
BRANCH="-r trilinos-release-10-0-branch"

# Source the 
#this is not the best thing to do, but we'll leave it for now
cd $HOME
source .bash_profile

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

echo
echo "Starting nightly Trilinos Release 10.0 testing on trilinos-test: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co $BRANCH -d $CHECKOUT_DIR Trilinos/cmake Trilinos/CTestConfig.cmake

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
