#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on geminga: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly

export TDD_DEBUG_VERBOSE=1
export TDD_FORCE_CMAKE_INSTALL=0
export TRIBITS_TDD_USE_SYSTEM_CTEST=1

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export TDD_HTTPS_PROXY="https://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"
export https_proxy="https://sonproxy.sandia.gov:80"

. ~/.bashrc

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
export CTEST_CONFIGURATION="default"
module load openmpi/1.10.0
module load gcc/5.2.0
module load valgrind/3.10.1

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`

echo "Configuration = $CONFIGURATION"
env

export OMP_NUM_THREADS=2

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

module unload valgrind
module unload gcc
module unload openmpi
# ===========================================================================
export CTEST_CONFIGURATION="nvcc_wrapper"
module load openmpi/1.10.0
module load gcc/4.9.2
module load cuda/7.5-gcc
module load nvcc-wrapper/gcc

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`

echo "Configuration = $CONFIGURATION"
env

export CUDA_LAUNCH_BLOCKING=1

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

module unload nvcc-wrapper
module unload cuda
module unload gcc
module unload openmpi
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on geminga: `date`"
echo
