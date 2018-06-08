#!/bin/bash

#
# Simple Jenkins driver for Trilinos builds where the build env is defined by
# sourcing the script:
#
#   source <trlinosBaseDir>/cmake/std/${CTEST_BUILD_CONFIGURATION_NAME}_env.sh
#
# and the build confiugration itself is defiled by pulling in the file:
#
#   <trilinosBaseDir>/cmake/std/$ENV{CTEST_BUILD_CONFIGURATION_NAME}.cmake
#
# in the CMake configure.
#
# To test this locally, first set up:
#
#   cd <some-base_dir>/
#   ln -s <trlinos-base-dir>/Trilinos .
#   mkdir SRC_AND_BUILD
#   cd SRC_AND_BUILD/
#   ln -s <trlinos-base-dir>/Trilinos .
#   cc ..
#
# then run:
#
#   env \
#     CTEST_BUILD_CONFIGURATION_NAME=<build-confg-name> \
#     Trilinos_PACKAGES=Teuchos \
#     CTEST_TEST_TYPE=Experimental \
#     CTEST_DO_SUBMIT=OFF \
#     CTEST_DO_UPDATES=OFF \
#     CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE \
#   <trlinos-base-dir>/Trilinos/cmake/ctest/drivers/sems_ci/ctest_std_driver.sh \
#     &> console.out
#
# To test the actual inner clone and update of Trilinos, don't do the inner
# symlink SRC_AND_BUILD/Trilinos and run with CTEST_DO_UPDATES=ON.
#


DRIVER_SCRIPT_DIR=`echo $0 | sed "s/\(.*\)\/.*\.sh/\1/g"`
echo "DRIVER_SCRIPT_DIR = '$DRIVER_SCRIPT_DIR'"

TRILINOS_DIR=`readlink -f ${DRIVER_SCRIPT_DIR}/../../../..`
echo "TRILINOS_DIR='${TRILINOS_DIR}'"

echo
echo "Loading env ${CTEST_BUILD_CONFIGURATION_NAME}_env.sh ..."
echo

source /etc/bashrc
source $TRILINOS_DIR/cmake/std/${CTEST_BUILD_CONFIGURATION_NAME}_env.sh

module list
echo
which cmake
echo
cmake --version
echo

export SUBDIR=SRC_AND_BUILD
if [ ! -e $SUBDIR ] ; then
  echo "Making $SUBDIR"
  mkdir $SUBDIR
fi

cd $SUBDIR/
echo "Current dir: $PWD"

export CTEST_DASHBOARD_ROOT=$PWD

env http_proxy= \
ctest -V -S $DRIVER_SCRIPT_DIR/ctest_std_driver.cmake
