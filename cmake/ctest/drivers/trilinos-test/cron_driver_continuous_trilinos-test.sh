#!/bin/bash

CTEST_EXE=/home/trilinos/cmake-2.8-rc4/bin/ctest
BASEDIR=/home/jmwille/TrilinosTestHarness/TrilinosDevelopment
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test

# Source the 
cd $HOME
source .bash_profile

if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

echo "_DAYOFWEEK=$_DAYOFWEEK"

echo
echo "Starting Trilinos continuous integration testing on trilinos-test: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake


echo
echo "Doing mpi debug continuous builds: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_continuous_mpi_debug_shared_trilinos-test.cmake -VV \
  &> $BASEDIR/ctest_linux_continuous_mpi_debug_shared_trilinos-test.out

echo
echo "Ending Trilinos continuous integration testing on trilinos-test. `date`"
echo


#/home/rabartl/mailmsg.py "Finished nightly Trilinos tests trilinos-test. http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
