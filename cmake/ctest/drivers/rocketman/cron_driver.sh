#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on rocketman: `date`"
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
module load sems-env
module load sems-cmake/3.5.2
module load sems-gcc/5.3.0
module load sems-openmpi/1.10.1
module load sems-superlu/4.3/base

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`

echo "Configuration = $CTEST_CONFIGURATION"
env

export OMP_NUM_THREADS=2

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

module unload sems-superlu/4.3/base
module unload sems-openmpi/1.10.1
module unload sems-gcc/5.3.0
module unload sems-cmake/3.5.2
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on rocketman: `date`"
echo
