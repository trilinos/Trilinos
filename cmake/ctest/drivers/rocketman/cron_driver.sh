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
export TDD_CTEST_TEST_TYPE=Nightly

# enable this to avoid clobbering any local changes you're making
#export TDD_IN_TESTING_MODE=ON

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

# update the Trilinos source
pushd $TRILINOS_SOURCE
./cmake/tribits/python_utils/gitdist --dist-no-color pull
popd

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
if [ -n "$DO_TPETRA_TESTING" ]; then
    ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_tpetra_rocketman.cmake
    ctest -S $SCRIPT_DIR/mpi_release_tpetra_deprecated_code_off_downstream_enabled_no_epetra.cmake
    ctest -S $SCRIPT_DIR/mpi_release_tpetra_deprecated_code_off_downstream_enabled.cmake
    ctest -S $SCRIPT_DIR/mpi_release_tpetra_deprecated_code_off_downstream_enabled_GO_int.cmake
    ctest -S $SCRIPT_DIR/ctest_linux_experimental_mpi_release_tpetra_performance_rocketman.cmake
else
    ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_rocketman.cmake
    ctest -S $SCRIPT_DIR/ctest_linux_experimental_mpi_release_avatar_rocketman.cmake
fi
popd

module unload muelu-gcc
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on $HOSTNAME: `date`"
echo
