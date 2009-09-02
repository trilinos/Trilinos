#!/bin/bash

CTEST_EXE=/home/jmwille/install/cmake/bin/ctest
BASEDIR=/space/jmwille/TrilinosTestHarness/TrilinosDevelopment
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/exetazo

# Source the 
cd $HOME
source .bash_profile

if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

echo "_DAYOFWEEK=$_DAYOFWEEK"

echo
echo "Starting Trilinos continuous integration testing on exetazo: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake


echo
echo "Doing mpi debug continuous builds: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_continuous_mpi_debug_shared_exetazo.cmake -VV \
  &> $BASEDIR/ctest_linux_continuous_mpi_debug_shared_exetazo.out

echo
echo "Ending Trilinos continuous integration testing on exetazo: `date`"
echo


#/home/rabartl/mailmsg.py "Finished nightly Trilinos tests exetazo: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
