#!/bin/bash

# OLCF Frontier MI250X

module purge
module load DefApps
module load PrgEnv-amd
module load amd/6.4.1
module load rocm/6.4.1
module swap cray-mpich cray-mpich/8.1.30
module load craype-accel-amd-gfx90a
module load cray-libsci
module load zlib
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load cray-parallel-netcdf
module load ninja
module load xpmem

export LD_LIBRARY_PATH="${CRAY_LD_LIBRARY_PATH}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

export HIP_ENABLE_DEFERRED_LOADING=0
export MPICH_GPU_SUPPORT_ENABLED=1
export TPETRA_ASSUME_GPU_AWARE_MPI=1
