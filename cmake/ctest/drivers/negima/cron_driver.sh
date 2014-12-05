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

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"

export TDD_FORCE_CMAKE_INSTALL=1

export CPATH=~/local/include:$CPATH
export LD_LIBRARY_PATH=~/local/lib:~/local/lib64:/usr/lib64:/lib64:$LD_LIBRARY_PATH
export PATH=~/bin:/opt/bin:~/local/bin:$PATH
export PKG_CONFIG_PATH=~/local:$PKG_CONFIG_PATH

module use-append ~/.modules
module load common_base
module load openmpi/1.8.3
module load gcc/4.9.2

export OMPI_MPICC=/home/aprokop/local/opt/gcc-4.9.2/bin/gcc
export OMPI_MPICXX=/home/aprokop/local/opt/gcc-4.9.2/bin/g++

env

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on negima: `date`"
echo
