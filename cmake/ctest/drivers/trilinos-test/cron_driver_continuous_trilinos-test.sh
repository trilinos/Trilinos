#!/bin/bash

CTEST_EXE=/home/trilinos/cmake/bin/ctest
EG_EXE=/home/trilinos/git/bin/eg
BASEDIR=/home/trilinos/dashboards/dev
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test
TRILINOS_REPOSITORY_LOCATION="software.sandia.gov:/space/git/Trilinos"

# Source the 
cd $HOME
source .bash_profile

if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

DATE=`date "+%Y%m%d"`

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
  echo Cloning the repository because none exists yet
  $EG_EXE clone $TRILINOS_REPOSITORY_LOCATION
fi


echo
echo "Doing continuous mpi optimized development shared library build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_continuous_mpi_development_opt_shared_trilinos-test.cmake -VV &> "CONTINUOUS_MPI_OPT_DEV_SHARED_$DATE.log"


echo
echo "Ending Trilinos continuous integration testing on trilinos-test. `date`"
echo
