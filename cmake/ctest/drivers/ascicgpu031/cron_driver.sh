#!/bin/bash --login

echo
echo "Starting nightly Trilinos development testing on `hostname`: `date`"
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


# Machine specific environment
#
export http_proxy=http://user:pass@proxy.sandia.gov:80/
export ftp_proxy=http://user:pass@proxy.sandia.gov:80/
export https_proxy=http://user:pass@proxy.sandia.gov:80/
export no_proxy=localhost,localnets,.sandia.gov
export TDD_HTTP_PROXY=$http_proxy
export TDD_HTTPS_PROXY=$https_proxy

. ~/.bashrc

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
export CTEST_CONFIGURATION="nvcc_wrapper"
module load sierra-devel/nvidia
module load sems-env
module unload sems-cmake
module load sems-cmake/3.17.1

# See Trilinos github issue #2115.
export OMPI_CXX=/home/csiefer/Trilinos/ascicgpu-testing/Trilinos/packages/kokkos/bin/nvcc_wrapper

# Remove colors (-fdiagnostics-color) from OMPI flags
# It may result in non-XML characters on the Dashboard
export OMPI_CFLAGS=`echo $OMPI_CFLAGS | sed 's/-fdiagnostics-color//'`

echo "OMPI_CXXFLAGS before $OMPI_CXXFLAGS"
export OMPI_CXXFLAGS=`echo $OMPI_CXXFLAGS | sed 's/-fdiagnostics-color//'`
echo "OMPI_CXXFLAGS after $OMPI_CXXFLAGS"
export CUDA_LAUNCH_BLOCKING=0
echo "Using cmake = `which cmake`, ctest = `which ctest`"

echo "Configuration = $CTEST_CONFIGURATION"
env

# Set the TMPDIR
mkdir -p /tmp/$USER
export TMPDIR=/tmp/$USER

# Machine independent cron_driver:
SCRIPT_DIR=`pwd`
$SCRIPT_DIR/../cron_driver.py

# Blitz the TMPDIR
rm -rf /tmp/$USER

module unload sems-cmake/3.17.1
module unload sems-env
module unload sierra-devel/nvidia
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on `hostname`: `date`"
echo
