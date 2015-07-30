#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on enigma: `date`"
echo

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

export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"
export TDD_FORCE_CMAKE_INSTALL=0
export CUDA_LAUNCH_BLOCKING=1
export OMP_NUM_THREADS=2
#export PATH="$PATH:/opt/intel/composer_xe_2015.1.133/bin/intel64"
#export LD_LIBRARY_PATH="/opt/intel/composer_xe_2015.1.133/compiler/lib/intel64:$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH="/usr/local/gcc/4.8.3/lib64:$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH="/opt/mpc/1.0.1/lib:$LD_LIBRARY_PATH"
#export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib:/usr/local/lib64"
# Machine independent cron_driver:
#
#openmpi-1.7-cuda6

. ~/.bashrc

# If you update the list of modules, go to ~/code/trilinos-test/trilinos/ and
# do "git pull". Otherwise, the tests could fail on the first night, as we
# would first run old cron_driver.sh and only then pull

# load OpenMPI 1.6.4 (standard in RHEL 7)
module load mpi/openmpi-x86_64

env


SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on enigma: `date`"
echo
