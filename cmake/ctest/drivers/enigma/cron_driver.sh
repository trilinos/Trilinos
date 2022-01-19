#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on $HOSTNAME: `date`"
echo

# Undefine the next line while making/testing local driver changes.  Otherwise, the nightly
# testing system will pull a fresh version of Trilinos and wipe out your changes.
# export TDD_IN_TESTING_MODE=1

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly
#export TDD_CTEST_TEST_TYPE=Experimental

export TDD_DEBUG_VERBOSE=1
export TRIBITS_TDD_USE_SYSTEM_CTEST=1

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
. /etc/profile

export TDD_HTTP_PROXY=$http_proxy
export TDD_HTTPS_PROXY=$https_proxy

export TDD_FORCE_CMAKE_INSTALL=0

. ~/.bashrc


# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
export MODULEPATH=$SCRIPT_DIR:$MODULEPATH
# Trilinos source repo
export TRILINOS_SOURCE=$SCRIPT_DIR/../../../..

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

module load muelu-gcc
module list
env

pushd $TRILINOS_SOURCE
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_debug_muelu_basker_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_debug_muelu_klu2_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_debug_muelu_extratypes_ei_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_serial_debug_muelu_extratypes_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_serial_release_muelu_experimental_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_mpi_release_no_serial_openmp_complex_experimental_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_mpi_release_openmp_no_epetra_no_int_complex_experimental_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_mpi_release_muelu_no_int_no_serial_openmp_experimental_enigma.cmake
ctest -S $SCRIPT_DIR/ctest_linux_nightly_mpi_release_muelu_no_int_openmp_experimental_enigma.cmake
popd

echo
echo "Ending nightly Trilinos development testing on $HOSTNAME: `date`"
echo
