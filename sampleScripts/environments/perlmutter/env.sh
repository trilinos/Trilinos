#!/bin/bash

# NERSC Perlmutter

module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load cray-parallel-netcdf
module load cray-libsci
module load cmake

# Because Intrepid2 won't build with GNU 13.2
module switch gcc-native/12.3

export HDF5_USE_FILE_LOCKING=FALSE

export CUDATOOLKIT_VERSION_STRING=${CRAY_CUDATOOLKIT_VERSION#*_}
export MPICH_GPU_SUPPORT_ENABLED=1
export TPETRA_ASSUME_GPU_AWARE_MPI=1

export CRAY_CPU_TARGET=x86-milan
export CRAY_ACCEL_TARGET=nvidia80

env_script_path=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export NVCC_WRAPPER=${env_script_path}/nvcc_wrapper
