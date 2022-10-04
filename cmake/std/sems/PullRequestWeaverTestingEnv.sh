# This script can be used to load the appropriate environment for the
# PR build on weaver using CUDA.

# usage: $ source PullRequestCUDA10.2.2TestingEnv.sh

#module load devpack/20210226/openmpi/4.0.5/gcc/7.2.0/cuda/10.2.2
#module load devpack/20190814/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
module load devpack/20210226/openmpi/4.0.5/gcc/7.2.0/cuda/10.2.2
module unload cmake
module load cmake/3.19.3
#current_dir=`dirname $BASH_SOURCE`
#Trilinos_dir=`realpath $current_dir/../../..`
#export OMPI_CXX=${Trilinos_dir}/packages/kokkos/bin/nvcc_wrapper
#echo ${OMPI_CXX}
#export OMPI_CC=`which gcc`
#export OMPI_FC=`which gfortran`
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
