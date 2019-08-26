#!/bin/csh

echo
echo "Starting nightly Trilinos development testing on dorksaber: `date`"
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
module load sems-env
module load sems-cmake/3.10.3
module load sems-gcc/5.3.0
module load sems-openmpi/1.10.1
module load sems-superlu/4.3/base

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
#setenv OMPI_CFLAGS="`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`"
#setenv OMPI_CXXFLAGS="`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`"

echo "Configuration = $CTEST_CONFIGURATION"
env

setenv OMP_NUM_THREADS 2



# Machine independent cron_driver:
setenv SCRIPT_DIR `dirname "$0"`
echo "SCRIPT_DIR = " $SCRIPT_DIR
$SCRIPT_DIR/../cron_driver.py

module unload sems-superlu/4.3/base
module unload sems-openmpi/1.10.1
module unload sems-gcc/5.3.0
module unload sems-cmake/3.10.3
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on dorksaber: `date`"
echo
