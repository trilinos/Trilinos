#
# This file is sourced in order to configure and build Trilinos with CUDA on the
# ATTB machines.
#
# NOTE: Set the env var TRILINOS_BASE_DIR to your location for Trilinos if it is not
# under $HOME/workspace.
#

if [ "$TRILINOS_BASE_DIR" == "" ] ; then
  export TRILINOS_BASE_DIR=$HOME/workspace
fi
#!/bin/bash
module purge
module load seacas/serial/20160328

export USE_CUDA=ON
NODE_TYPE=CUDA
export NODE_TYPE

echo using compiler stack $COMPILER to build $BUILD_TYPE code

module load devpack/openmpi/1.10.0/gcc/4.8.4/cuda/7.5.18
export OMPI_CXX=$TRILINOS_BASE_DIR/Trilinos/packages/kokkos/config/nvcc_wrapper 
export USE_CUDA=ON
module load blas/gcc/4.8.4
module load lapack/3.5.0/gcc/4.8.4
export LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-L${BLAS_ROOT}/lib;-lblas;-lgfortran"


module swap yaml-cpp/0.3.0 yaml-cpp/20170104 
if [ $? ]; then module load  yaml-cpp/20170104; fi
export COMPILER
unset ATTB_ENV

export NETCDF_BASE_DIR=$NETCDF_ROOT
export HDF_BASE_DIR=$HDF5_ROOT
export BOOST_BASE_DIR=$BOOST_ROOT
export BUILD_TYPE=DEBUG

echo BUILD_TYPE = ${BUILD_TYPE}
echo Using openmp = ${USE_OPENMP}
