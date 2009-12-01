#!/bin/bash

# Set what ctest to use
CTEST_EXE=/usr/local/bin/ctest


echo
echo "Starting nightly Trilinos testing on gabriel: `date`"
echo

#BASEDIR=`echo $0 | sed "s/\(.*\)\/Trilinos\/cmake\/ctest\/.*\.sh/\1/g"`
BASEDIR=/mnt/disk2/rabartl/CDash/Trilinos.base
echo "BASEDIR = $BASEDIR"
DRIVER_SCRIPT_DIR=$BASEDIR/Trilinos/cmake/ctest/drivers/gabriel
echo "DRIVER_SCRIPT_DIR = $DRIVER_SCRIPT_DIR"


echo
echo "Setting up thet environment ..."
echo

cd $DRIVER_SCRIPT_DIR
source bash_profile
echo $PATH
which eg


echo
echo "Updating Trilinos driver code: `date`"
echo

cd $BASEDIR/Trilinos
eg pull

cd $BASEDIR
  
echo
echo "Doing mpi debug build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_mpi_debug_gabriel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_mpi_debug_gabriel.out


echo
echo "Kill any remaining 'orted' processes: `date`"
echo

killall -s 9 orted


echo
echo "Doing serial release build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_release_gabriel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_release_gabriel.out


echo
echo "Doing serial debug boost tracing build: `date`"
echo

time ${CTEST_EXE} -S $DRIVER_SCRIPT_DIR/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.cmake -VV \
  &> $BASEDIR/ctest_linux_nightly_serial_debug_boost_tracing_gabriel.out


echo
echo "Ending nightly Trilinos testing on gabriel: `date`"
echo


/home/rabartl/mailmsg.py "Finished nightly Trilinos CMake tests gabriel: http://trilinos-dev.sandia.gov/cdash/index.php?project=Trilinos"
