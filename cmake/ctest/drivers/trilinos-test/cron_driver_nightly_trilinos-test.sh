#!/bin/bash

CTEST_EXE=/home/trilinos/cmake-2.8-rc4/bin/ctest
BASEDIR=/home/bmpersc/nightly/Trilinos.base/release_10
BASEDATADIR=/home/bmpersc/nightly
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test
BRANCH="-r trilinos-release-10-0-branch"

export CMAKE_LIBRARY_PATH="/home/trilinos/tpl/gcc4.1.2/exodusII_4.84/lib:/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/lib"
export CMAKE_INCLUDE_PATH="/home/trilinos/tpl/gcc4.1.2/exodusII_4.84/include:/home/trilinos/tpl/gcc4.1.2/netcdf_4.0/include"
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
cvs -q -d :ext:software:/space/CVS co $BRANCH Trilinos/cmake Trilinos/CTestConfig.cmake

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
