#!/bin/csh

echo
echo "Starting nightly Trilinos development testing on lightsaber: `date`"
echo

#
# TrilinosDriver settings:
#

setenv TDD_PARALLEL_LEVEL 2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
setenv TDD_CTEST_TEST_TYPE Nightly


# Machine specific environment
#

setenv TDD_HTTP_PROXY "http://wwwproxy.sandia.gov:80"
setenv TDD_HTTPS_PROXY "https://wwwproxy.sandia.gov:80"
setenv http_proxy "http://wwwproxy.sandia.gov:80"
setenv https_proxy "https://wwwproxy.sandia.gov:80"
setenv TDD_FORCE_CMAKE_INSTALL 1
setenv TDD_DEBUG_VERBOSE 1

source ~/.cshrc


# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
setenv CTEST_CONFIGURATION "default"
module purge
module load sems-gcc/10.1.0
module load sems-openmpi/4.0.5
module load sems-cmake/3.21.1
module load sems-superlu/4.3
module load sems-zlib/1.2.11
module load sems-boost/1.74.0
module load sems-netcdf-c/4.7.3

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
#setenv OMPI_CFLAGS="`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`"
#setenv OMPI_CXXFLAGS="`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`"

echo "Configuration = $CTEST_CONFIGURATION"
env

setenv OMP_NUM_THREADS 2

# Update Avatar 
(cd /home/nightlyTesting/avatar; git pull --rebase )

# Set variables to work aroun TriBITS problems
#setenv TDD_FORCE_CMAKE_INSTALL 0
setenv TRIBITS_TDD_USE_SYSTEM_CTEST 1

# Machine independent cron_driver:
setenv SCRIPT_DIR `dirname "$0"`
echo "SCRIPT_DIR = " $SCRIPT_DIR
$SCRIPT_DIR/../cron_driver.py

module load sems-netcdf-c
module load sems-boost
module load sems-zlib
module unload sems-superlu
module unload sems-cmake
module unload sems-openmpi
module unload sems-gcc
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on lightsaber: `date`"
echo
