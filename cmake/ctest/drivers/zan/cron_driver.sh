#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on zan: `date`"
echo

# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=1

# Trilinos settings:
#

#export CTEST_TEST_TYPE=Experimental

#export CTEST_DO_SUBMIT=FALSE

#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

export PATH=/usr/kerberos/bin:/opt/intel/11.0.074/bin/intel64:/opt/intel/11.0.074/bin/intel64:/usr/local/bin:/bin:/usr/bin:/usr/local/cuda/bin:/home/jmwille/bin:$PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/11.0.074/ipp/em64t/sharedlib:/opt/intel/11.0.074/mkl/lib/em64t:/opt/intel/11.0.074/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib:/opt/intel/11.0.074/lib/intel64:/opt/intel/11.0.074/ipp/em64t/sharedlib:/opt/intel/11.0.074/mkl/lib/em64t:/opt/intel/11.0.074/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib:/opt/intel/11.0.074/lib/intel64

#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on zan: `date`"
echo
