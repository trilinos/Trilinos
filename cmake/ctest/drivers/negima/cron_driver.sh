#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on negima: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=Nightly

export TDD_FORCE_CMAKE_INSTALL=1

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"

. ~/.bashrc

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# ===========================================================================
export CTEST_CONFIGURATION="default"
module load openmpi/1.10.0
module load gcc/5.2.0
module load valgrind/3.10.1

echo "Configuration = $CONFIGURATION"
env

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

module unload valgrind
module unload gcc
module unload openmpi
# ===========================================================================
export CTEST_CONFIGURATION="nvcc_wrapper"
module load openmpi/1.10.0
module load gcc/4.8.4
module load cuda/6.5-gcc
module load nvcc-wrapper/gcc

echo "Configuration = $CONFIGURATION"
env

# Machine independent cron_driver:
SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

module unload nvcc-wrapper
module unload cuda
module unload gcc
module unload openmpi
# ===========================================================================

echo
echo "Ending nightly Trilinos development testing on negima: `date`"
echo
