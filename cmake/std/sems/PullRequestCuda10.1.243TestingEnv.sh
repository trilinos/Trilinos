# This script can be used to load the appropriate environment for the
# PR build on vortex using CUDA.

# usage: $ source PullRequestCUDA10.1.105TestingEnv.sh

module load StdEnv
module load sparc-dev/cuda-10.1.243_gcc-7.3.1_spmpi-rolling
#export OMPI_CXX=`which g++`
current_dir=`dirname $BASH_SOURCE`
Trilinos_dir=`realpath $current_dir/../../..`
export OMPI_CXX=${Trilinos_dir}/packages/kokkos/bin/nvcc_wrapper
export OMPI_CC=`which gcc`
export OMPI_FC=`which gfortran`
export CUDA_LAUNCH_BLOCKING=1
export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1

# Use manually installed cmake and ninja to try to avoid module loading
# problems (see TRIL-208)
export PATH=/projects/atdm_devops/vortex/cmake-3.17.2/bin:$PATH
export PATH=/projects/atdm_devops/vortex/ninja-fortran-1.8.2:$PATH

