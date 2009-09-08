#!/bin/bash

CTEST_EXE=/home/trilinos/cmake/bin/ctest
BASEDIR=/home/bmpersc/nightly/Trilinos.base
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test

# Source the 
#this is not the best thing to do, but we'll leave it for now
cd $HOME
source .bash_profile

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

echo
echo "Starting nightly Trilinos testing on trilinos-test: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

echo
echo "Doing mpi optimized release 10.0 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_release_trilinos-test.cmake -VV &> "MPI_RELEASE_10.0_$DATE.log"

echo
echo "Doing mpi optimized release 10.0  shared library build: `date`"
echo "Temporarily disabled do to hang on submit"

#time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_release_shared_trilinos-test.cmake -VV &> "MPI_RELEASE_10.0_SHARED_$DATE.log"

echo
echo "Doing serial release 10.0 build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_release_trilinos-test.cmake -VV &> "SERIAL_RELEASE_10.0_$DATE.log"

#These are development runs. Eventually they will need to be enabled.

#echo
#echo "Doing mpi optimized development build: `date`"
#echo
#
#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/drivers/s909348/ctest_mac_nightly_mpi_opt_s909348.cmake -VV &> "MPI_OPT_DEV_$DATE.log"
#
#echo
#echo "Doing serial debug development build: `date`"
#echo
#
#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/drivers/s909348/ctest_mac_nightly_serial_debug_s909348.cmake -VV &> "SERIAL_DEBUG_DEV_$DATE.log"

echo
echo "Ending nightly Trilinos testing on trilinos-test: `date`"
echo
