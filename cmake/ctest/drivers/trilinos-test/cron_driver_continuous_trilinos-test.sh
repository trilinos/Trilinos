#!/bin/bash

CTEST_EXE=/home/trilinos/cmake/bin/ctest
EG_EXE=/home/trilinos/git/bin/eg
BASEDIR=/home/jmwille/TrilinosTestHarness/TrilinosDevelopment
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

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
echo "Doing mpi debug continuous builds: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_continuous_mpi_debug_shared_trilinos-test.cmake -VV \
  &> $BASEDIR/ctest_linux_continuous_mpi_debug_shared_trilinos-test.out

echo
echo "Ending Trilinos continuous integration testing on trilinos-test. `date`"
echo


#/home/rabartl/mailmsg.py "Finished nightly Trilinos tests trilinos-test. http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
