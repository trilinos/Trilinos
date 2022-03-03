#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on $HOSTNAME: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
#export TDD_CTEST_TEST_TYPE=Specialized
export TDD_CTEST_TEST_TYPE=Nightly

export TDD_DEBUG_VERBOSE=1
export TDD_FORCE_CMAKE_INSTALL=0
export TRIBITS_TDD_USE_SYSTEM_CTEST=1

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
. /etc/profile

export TDD_HTTP_PROXY=$http_proxy
export TDD_HTTPS_PROXY=$https_proxy


. ~/.bashrc

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
export MODULEPATH=$SCRIPT_DIR:$MODULEPATH
# Trilinos source repo
export TRILINOS_SOURCE=$SCRIPT_DIR/../../../..

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
export CTEST_CONFIGURATION="default"
module load muelu-gcc
module list

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`

echo "Configuration = $CTEST_CONFIGURATION"
env

pushd $TRILINOS_SOURCE
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_trappist.gcc.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_no_int_no_serial_openmp_trappist.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_experimental_openmp_refactor_muelu_trappist.gcc.cmake
popd

module unload muelu-gcc
# ===========================================================================
export CTEST_CONFIGURATION="clang"
module load muelu-clang
module list

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`

echo "Configuration = $CTEST_CONFIGURATION"
env

pushd $TRILINOS_SOURCE
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_trappist.clang.cmake
popd

module unload muelu-clang
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on $HOSTNAME: `date`"
echo
