#!/bin/bash

source /opt/casldev/env/casl_dev_env.sh

BASEDIR=/home/ogb/Dashboards
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/u233
TRILINOS_REPOSITORY_LOCATION="cgbaker@software.sandia.gov:/space/git/Trilinos"

export TDD_PARALLEL_LEVEL=2
export TDD_CTEST_TEST_TYPE=Experimental

export CTEST_DROP_SITE=casl-dev.ornl.gov
export CTEST_DROP_LOCATION="/CDash/submit.php?project=Trilinos"

if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

echo "_DAYOFWEEK=$_DAYOFWEEK"

if [ "$_DAYOFWEEK" == "Saturday" ] ; then
  _RUN_REGULAR_TESTS=1
elif [ "$_DAYOFWEEK" == "Sunday" ] ; then
  _RUN_REGULAR_TESTS=1 
else
  _RUN_REGULAR_TESTS=1
fi

echo "_RUN_REGULAR_TESTS=$_RUN_REGULAR_TESTS"

#exit

echo
echo "Starting nightly Trilinos testing on u233: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

cd $BASEDIR
if [ -d Trilinos ]; then
  echo Doing an update of existing directory
  cd Trilinos
  eg pull
  cd ..
else
  echo Cloning the repository because none exists yets
  eg clone $TRILINOS_REPOSITORY_LOCATION
fi


#
# Weekday tests
#

if [ "$_RUN_REGULAR_TESTS" == "1" ] ; then

echo
echo "Running regular tests ..."
echo

echo
echo "Doing serial debug intel build: `date`"
echo

time ctest -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_icpc_u233.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_icpc_u233.out

echo
echo "Doing serial release intel build: `date`"
echo

time ctest -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_release_icpc_u233.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_release_icpc_u233.out

fi

echo
echo "Ending nightly Trilinos testing on u233: `date`"
echo

echo "Finished nightly Trilinos CMake tests u233: http://casl-dev.ornl.gov/CDash/index.php?project=Trilinos" | mailx -s "Nightly CTest: billmp1" bakercg@ornl.gov
#echo "Finished nightly Trilinos CMake tests u233: http://casl-dev.ornl.gov/CDash/index.php?project=Trilinos" | mailx -s "Nightly CTest: billmp1" -c rabartl@sandia.gov bakercg@ornl.gov
