# This script can be used to load the appropriate environment for the
# PR build on ride using CUDA.

# usage: $ source PullRequestCUDA9.2TestingEnv.sh

#No SEMS NFS mount on ride
#source /projects/sems/modulefiles/utils/sems-modules-init.sh
module load git/2.10.1
module load devpack/20180521/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88
module swap openblas/0.2.20/gcc/7.2.0 netlib/3.8.0/gcc/7.2.0
export OMPI_CXX=`which g++`
export OMPI_CC=`which gcc`
export OMPI_FC=`which gfortran`
export CUDA_LAUNCH_BLOCKING=1
export NVCC_WRAPPER_DEFAULT_COMPILER=$(WORKSPACE)/Trilinos/packages/kokkos/bin/nvcc_wrapper

# Use manually installed cmake and ninja to try to avoid module loading
# problems (see TRIL-208)
export PATH=/ascldap/users/rabartl/install/white-ride/cmake-3.11.2/bin:/ascldap/users/rabartl/install/white-ride/ninja-1.8.2/bin:$PATH
