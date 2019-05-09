#!/bin/bash

echo
echo "Starting nightly Trilinos development testing on typhon: `date`"
echo

#
# TrilinosDriver settings:
#

export TDD_PARALLEL_LEVEL=1

# Trilinos settings:
#

# Submission mode for the *TrilinosDriver* dashboard
export TDD_CTEST_TEST_TYPE=${JENKINS_JOB_TYPE}

#export CTEST_DO_SUBMIT=FALSE
#export CTEST_START_WITH_EMPTY_BINARY_DIRECTORY=FALSE

# Machine specific environment
#
module purge
module load sems-env
module load kokkos-env

module load sems-python/2.7.9
module load sems-cmake/3.10.3
module load sems-git/2.1.3
module load sems-${COMPILER_SUFFIX}

export TRIBITS_TDD_USE_SYSTEM_CTEST=1
export TDD_HTTP_PROXY="http://sonproxy.sandia.gov:80"
export http_proxy="http://sonproxy.sandia.gov:80"
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
export OMP_NUM_THREADS=2

# Machine independent cron_driver:
#
#openmpi-1.7-cuda6

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`

if [ "${JENKINS_DO_CUDA}" == 'ON' ]; then
  module load kokkos-${CUDA_SUFFIX}
  module load kokkos-${MPI_SUFFIX}/cuda
  export OMPI_CXX=$WORKSPACE/Trilinos/packages/kokkos/config/nvcc_wrapper
else
  module load sems-${MPI_SUFFIX}
fi

module load sems-${BOOST_SUFFIX}/base
module load sems-${HDF5_SUFFIX}/parallel
module load sems-${NETCDF_SUFFIX}/exo_parallel
module load sems-${ZLIB_SUFFIX}/base

module list

$SCRIPT_DIR/../cron_driver.py

echo
echo "Ending nightly Trilinos development testing on typhon: `date`"
echo
