#!/bin/bash

CTEST_EXE=/home/trilinos/bin/ctest
EG_EXE=/home/trilinos/bin/eg
BASEDIR=/home/trilinos/dashboards/cont-dev
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/trilinos-test2
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
echo "Starting Trilinos continuous integration testing on trilinos-test2: `date`"
echo


export LD_LIBRARY_PATH="/home/trilinos/gcc4.7.2/base/lib64:$LD_LIBRARY_PATH"

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

time ${CTEST_EXE} -j 14 -S $DRIVER_SCRIPT_DIR/ctest_linux_continuous_mpi_opt_shared_trilinos-test2.cmake &> "CONTINUOUS_MPI_OPT_DEV_SHARED_$DATE.log"


echo
echo "Ending Trilinos continuous integration testing on trilinos-test2. `date`"
echo

