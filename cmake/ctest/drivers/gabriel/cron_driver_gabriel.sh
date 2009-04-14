#!/bin/bash

# Source the 
cd $HOME
source .bash_profile

CTEST_EXE=/usr/local/bin/ctest

if [ "$_DAYOFWEEK" == "" ] ; then
  _DAYOFWEEK=`date +%A`
fi

echo "_DAYOFWEEK=$_DAYOFWEEK"

if [ "$_DAYOFWEEK" == "Saturday" ] ; then
  _RUN_REGULAR_TESTS=1
  _RUN_COVERAGE_TESTS=1
  _RUN_MEMCHECK_TESTS=0
elif [ "$_DAYOFWEEK" == "Sunday" ] ; then
  _RUN_REGULAR_TESTS=0 
  _RUN_COVERAGE_TESTS=0
  _RUN_MEMCHECK_TESTS=1
else
  _RUN_REGULAR_TESTS=1
  _RUN_COVERAGE_TESTS=0
  _RUN_MEMCHECK_TESTS=0
fi

echo "_RUN_REGULAR_TESTS=$_RUN_REGULAR_TESTS"
echo "_RUN_COVERAGE_TESTS=$_RUN_COVERAGE_TESTS"
echo "_RUN_MEMCHECK_TESTS=$_RUN_MEMCHECK_TESTS"

#exit


#export Trilinos_PACKAGES=Teuchos


echo
echo "Starting nightly Trilinos testing on gabriel: `date`"
echo


echo
echo "Checking out just the skeleton cmake/ctest code: `date`"
echo

BASEDIR=/home/rabartl/PROJECTS/dashboards/Trilinos.base

cd $BASEDIR
cvs -q -d :ext:software:/space/CVS co Trilinos/cmake Trilinos/CTestConfig.cmake


#
# Weekday tests
#

if [ "$_RUN_REGULAR_TESTS" == "1" ] ; then
  
  echo
  echo "Running weekday tests ..."
  echo
  
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

fi

echo
echo "Ending nightly Trilinos testing on gabriel: `date`"
echo


/home/rabartl/mailmsg.py "Finished nightly Trilinos tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
