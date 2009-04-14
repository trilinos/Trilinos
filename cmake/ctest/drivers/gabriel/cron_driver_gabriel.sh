#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

CTEST_EXE=/usr/local/bin/ctest


echo
echo "Starting nightly Trilinos testing on gabriel: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

BASEDIR=/home/rabartl/PROJECTS/dashboards/Trilinos.base

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake

  
echo
echo "Doing mpi optimized build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/drivers/gabriel/ctest_linux_nightly_mpi_optimized_gabriel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_optimized_gabriel.out


echo
echo "Doing serial debug build: `date`"
echo

time ${CTEST_EXE} -S $BASEDIR/Trilinos/cmake/ctest/drivers/gabriel/ctest_linux_nightly_serial_debug_gabriel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_gabriel.out


echo
echo "Ending nightly Trilinos testing on gabriel: `date`"
echo


/home/rabartl/mailmsg.py "Finished nightly Trilinos tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
