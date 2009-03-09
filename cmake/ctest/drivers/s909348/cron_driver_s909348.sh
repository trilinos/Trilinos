#!/bin/bash

# Source the 
#this is not the best thing to do, but we'll leave it for now
cd $HOME
source .bash_profile

#get the date for use in log files
DATE=`date "+%m-%d-%Y"`

CTEST_EXE=/Users/bmpersc/bin/cmake/bin/ctest

echo
echo "Starting nightly Trilinos testing on s909348: `date`"
echo

echo
echo "Checking out just the drivers: `date`"
echo

BASEDIR=/Users/bmpersc/nightly/Trilinos.base

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

#echo
#echo "Doing dependency checking-only build: `date`"
#echo
#
#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_package_deps_godel.cmake -VV

echo
echo "Doing serial performance build: `date`"
echo "temporarily disabled"
echo

#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_serial_performance_godel.cmake -VV

echo
echo "Doing mpi optimized build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_mac_nightly_mpi_release_s909348.cmake -VV &> "MPI_RELEASE_$DATE.log"

echo
echo "Doing serial debug build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_mac_nightly_serial_debug_s909348.cmake -VV &> "SERIAL_DEBUG_$DATE.log"

echo
echo "Doing mpi optimized shared library build: `date`"
echo "temporarily disabled"
echo

#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_optimized_shared_godel.cmake -VV

echo
echo "Doing mpi optimized zoltan c-only build: `date`"
echo "temporarily disabled"
echo

#time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/ctest_linux_nightly_mpi_opt_zoltan_c_godel.cmake -VV

echo
echo "Ending nightly Trilinos testing on s909348: `date`"
echo

#maybe I should replicate this. but later :)
#/home/rabartl/mailmsg.py "Finished nightly Trilinos tests godel: http://trilinos.sandia.gov/cdash/index.php?project=Trilinos"
