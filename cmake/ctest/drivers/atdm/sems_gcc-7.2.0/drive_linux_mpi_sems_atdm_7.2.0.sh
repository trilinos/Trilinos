#!/bin/bash

# This is a simple driver script that uses a cloned and updated Trilinos git
# repo to get the driver scripts which are run.  The ctest -S driver script
# will clone and update another Trilinos repo which will be used to drive the
# build.
#
# Just run this script from the base location where you want to create create
# the directories and that is all.

DRIVER_SCRIPT_DIR_REL=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
DRIVER_SCRIPT_DIR=`readlink -f ${DRIVER_SCRIPT_DIR_REL}`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../../..`
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

source /etc/bashrc
source $TRILINOS_DIR/cmake/std/sems/atdm/load_atdm_7.2_dev_env.sh

echo
module list

echo
echo "Some of the set env vars:"
set | grep "SEMS_.*_ROOT"
set | grep TRILINOS_SEMS_DEV_ENV_LOADED

TEST_DIR=SRC_AND_BUILD
if [ ! -e $TEST_DIR ] ; then
  echo "Creating $TEST_DIR"
  mkdir $TEST_DIR
fi
cd $TEST_DIR/

export Trilinos_REPOSITORY_LOCATION=https://github.com/trilinos/Trilinos.git

if [ "${Trilinos_CTEST_DO_ALL_AT_ONCE}" == "" ] ; then
  export Trilinos_CTEST_DO_ALL_AT_ONCE=TRUE
fi

env CTEST_DASHBOARD_ROOT=$PWD \
  ctest -V -S $DRIVER_SCRIPT_DIR/ctest_linux_mpi_sems_atdm_7.2.0.cmake
