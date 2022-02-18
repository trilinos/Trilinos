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

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`


# Trilinos source repo
export TRILINOS_SOURCE=$SCRIPT_DIR/../../../..

# folder with the machine specific build info
export BUILDS_DIR=$TRILINOS_SOURCE/cmake/ctest/drivers/$HOSTNAME

# So we can load the env module specific to the build
export MODULEPATH=$BUILDS_DIR:$MODULEPATH
echo "Loading modules from $MODULEPATH"

# update the Trilinos source
pushd $TRILINOS_SOURCE
./cmake/tribits/python_utils/gitdist --dist-no-color pull
popd


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

case "$(date +%a)" in
    Sun|Tue|Thu)
        heavy_memory_tests=true
        coverage=false
        ;;
    Mon|Wed)
        heavy_memory_tests=false
        coverage=true
        ;;
    *)
        heavy_memory_tests=false
        coverage=false
        ;;
esac

pushd $TRILINOS_SOURCE
ctest -S $BUILDS_DIR/ctest_linux_nightly_serial_debug_muelu_tpetra_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_serial_debug_muelu_epetra_geminga.cmake
# ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_complex_muelu_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_muelu_tpetra_no_int_no_serial_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_muelu_no_epetra_no_serial_openmp_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_tpetra_no_int_experimental_geminga.cmake

if [ "$heavy_memory_tests" = true ] ; then
    ctest -S $BUILDS_DIR/ctest_linux_nightly_serial_debug_valgrind_muelu_geminga.cmake
    ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_debug_muelu_geminga.cmake
else
    ctest -S $BUILDS_DIR/ctest_linux_nightly_serial_debug_muelu_geminga.cmake
fi

if [ "$coverage" = true ] ; then
   ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_debug_muelu_coverage_geminga.cmake
fi

popd

# ===========================================================================
export CTEST_CONFIGURATION="nvcc_wrapper"

# See Trilinos github issue #2115.
export OMPI_CXX=/home/jhu/code/trilinos-test/trilinos/packages/kokkos/bin/nvcc_wrapper

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`

echo "OMPI_CXXFLAGS before $OMPI_CXXFLAGS"
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`
echo "OMPI_CXXFLAGS after $OMPI_CXXFLAGS"

echo "Configuration = $CTEST_CONFIGURATION"
env

pushd $TRILINOS_SOURCE
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_muelu_kokkos_refactor_cuda_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_muelu_amgx_cuda_geminga.cmake
ctest -S $BUILDS_DIR/ctest_linux_nightly_mpi_release_muelu_cuda_no_uvm_geminga.cmake
popd

module unload muelu-gcc
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on $HOSTNAME: `date`"
echo
