#!/bin/bash

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
source /projects/modulefiles/utils/sems-modules-init.sh
source /projects/modulefiles/utils/kokkos-modules-init.sh

module load python/2.7.9
module load cmake/2.8.12
module load git/2.1.3

export TRIBITS_TDD_USE_SYSTEM_CTEST=1
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
export OMP_NUM_THREADS=2

# Machine independent cron_driver:
#
#openmpi-1.7-cuda6

SCRIPT_DIR=`cd "\`dirname \"$0\"\`";pwd`

if [ "${JENKINS_DO_CUDA}" == 'ON' ]; then
  module load ${CUDA_SUFFIX}
  module load ${COMPILER_SUFFIX}/${MPI_SUFFIX}/${CUDA_SUFFIX}
  export OMPI_CXX=$WORKSPACE/Trilinos/packages/kokkos/config/nvcc_wrapper
else
  module load $SEMS_MODULE_ROOT/rhel6-x86_64/sems/compiler/${COMPILER_SUFFIX}/${MPI_SUFFIX}
fi

module swap ${COMPILER_SUFFIX}/base ${COMPILER_SUFFIX}/base

module load ${BOOST_SUFFIX}/${COMPILER_SUFFIX}/base
module load ${HDF5_SUFFIX}/${COMPILER_SUFFIX}/${MPI_SUFFIX}
module load ${NETCDF_SUFFIX}/${COMPILER_SUFFIX}/parallel
module load ${ZLIB_SUFFIX}/${COMPILER_SUFFIX}/base

module list

$SCRIPT_DIR/../cron_driver.py

