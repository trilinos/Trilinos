# This script can be used to load the appropriate environment for the
# PR build on ride using CUDA.

# usage: $ source PullRequestCUDA10.1.105uvmOffTestingEnv.sh

module load git/2.10.1
module load devpack/20190404/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
module swap openblas/0.3.4/gcc/7.4.0 netlib/3.8.0/gcc/7.2.0
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
export PATH=/home/atdm-devops-admin/tools/ride/cmake-3.17.2/bin/:/ascldap/users/rabartl/install/white-ride/ninja-1.8.2/bin:$PATH
