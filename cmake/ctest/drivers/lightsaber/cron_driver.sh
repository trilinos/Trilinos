#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on lightsaber: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly


# Machine specific environment
#

export TDD_HTTP_PROXY="http://wwwproxy.sandia.gov:80"
export TDD_HTTPS_PROXY="https://wwwproxy.sandia.gov:80"
export http_proxy="http://wwwproxy.sandia.gov:80"
export https_proxy="https://wwwproxy.sandia.gov:80"
export TDD_FORCE_CMAKE_INSTALL=1
export TDD_DEBUG_VERBOSE=1

. /etc/profile
source ~/.bashrc


# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`

# Trilinos source repo
export TRILINOS_SOURCE=$SCRIPT_DIR/../../../..

# folder with the machine specific build info
export BUILDS_DIR=$TRILINOS_SOURCE/cmake/ctest/drivers/$HOSTNAME

# OneAPI
export MODULEPATH="$MODULEPATH":/opt/intel/oneapi/modulefiles:/opt/apps/modulefiles

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
# GCC family
echo "GREP: *** GCC Family Tests ***"
export CTEST_CONFIGURATION="default"
module purge
module load sems-gcc/10.1.0
module load sems-openmpi/4.0.5
module load sems-cmake
module load sems-superlu/4.3
module load sems-zlib
module load sems-boost
module load sems-hdf5
module load sems-netcdf-c
module load sems-parallel-netcdf

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
#setenv OMPI_CFLAGS="`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`"
#setenv OMPI_CXXFLAGS="`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`"

echo "Configuration = $CTEST_CONFIGURATION"
env

export OMP_NUM_THREADS=2

# Update Avatar
(cd /home/nightlyTesting/avatar; git pull --rebase )

# Set variables to work aroun TriBITS problems
#setenv TDD_FORCE_CMAKE_INSTALL 0
export TRIBITS_TDD_USE_SYSTEM_CTEST=1

# Actually run stuff
ctest -S $BUILDS_DIR/ctest_linux_experimental_mpi_release_avatar_lightsaber.cmake
ctest -S $BUILDS_DIR/ctest_linux_experimental_mpi_release_float_lightsaber.cmake

module unload sems-parallel-netcdf
module unload sems-netcdf-c
module unload sems-hdf5
module unload sems-boost
module unload sems-zlib
module unload sems-superlu
module unload sems-cmake
module unload sems-openmpi
module unload sems-gcc
# ===========================================================================


echo
echo "Ending nightly Trilinos development testing on lightsaber: `date`"
echo
