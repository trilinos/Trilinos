#!/bin/bash

#necessary because the testbeds don't setup modules by
#default on the login node or compute nodes.
source /etc/profile.d/modules.sh

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
module load python/2.7.12
export PATH=/home/rabartl/install/white-ride/cmake-3.11.2/bin:$PATH
module load git/2.10.1

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
  module load ${MPI_SUFFIX}/${COMPILER_SUFFIX}/${CUDA_SUFFIX}
  export OMPI_CXX=$WORKSPACE/Trilinos/packages/kokkos/config/nvcc_wrapper
else
  module load ${MPI_SUFFIX}/${COMPILER_SUFFIX}
fi

module load ${BOOST_SUFFIX}/${MPI_SUFFIX}/${COMPILER_SUFFIX}/${CUDA_SUFFIX}
module load ${HDF5_SUFFIX}/${MPI_SUFFIX}/${COMPILER_SUFFIX}/${CUDA_SUFFIX}
module load ${NETCDF_SUFFIX}/${MPI_SUFFIX}/${COMPILER_SUFFIX}/${CUDA_SUFFIX}
module load ${ZLIB_SUFFIX}
module load ${SUPERLU_SUFFIX}/${COMPILER_SUFFIX}
#module load ${BLAS_SUFFIX}/${COMPILER_SUFFIX}
#module load ${LAPACK_SUFFIX}/${COMPILER_SUFFIX}

module list

$SCRIPT_DIR/../cron_driver.py

