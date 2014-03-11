#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on s956330: `date`"
echo

# TrilinosDriver settings:
#

export TDD_GIT_EXE=/usr/bin/git
export TDD_PARALLEL_LEVEL=2

# Trilinos settings:
#

export TDD_CTEST_TEST_TYPE=Experimental

#export CTEST_DO_SUBMIT=FALSE

#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

#export Trilinos_PACKAGES=Teuchos

# Machine specific environment:
#

#export PYTHONPATH=/Users/jmwille/install/lib/python2.5/site-packages

export PATH=/projects/sems/utilities/bin:/projects/sems/compilers/clang/3.4/bin:/projects/sems/compilers/gcc/4.8.2/bin:/projects/sems/compilers/gcc/openmpi_1.6.5_gcc_4.8.2/bin:$PATH

export DYLD_LIBRARY_PATH=/projects/sems/compilers/clang/3.4/lib:/projects/sems/compilers/gcc/4.8.2/lib:/projects/sems/compilers/gcc/openmpi_1.6.5_gcc_4.8.2/lib:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=/projects/sems/compilers/clang/3.4/lib:/projects/sems/compilers/gcc/4.8.2/lib:/projects/sems/compilers/gcc/openmpi_1.6.5_gcc_4.8.2/lib:$LD_LIBRARY_PATH
export PYTHONPATH="/projects/sems/tpls/gcc_4.8.2/arch_x86/numpy_1.8.0/lib/python2.7/site-packages"

# Machine independent cron_driver:
#

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`
$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on s956330: `date`"
echo
