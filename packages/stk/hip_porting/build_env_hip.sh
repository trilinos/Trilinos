#!/bin/sh

module load craype-accel-amd-gfx908
module load PrgEnv-cray
module load rocm
module load boost
module load cray-hdf5
module load cray-netcdf
module load cray-parallel-netcdf

export MPIR_CVAR_GPU_EAGER_DEVICE_MEM=0
export MPICH_GPU_SUPPORT_ENABLED=1
export PATH=$HOME/Trilinos/builds:$PATH

cd ..
git checkout -- packages
git apply builds/hip_kokkos.patch
git apply builds/hip_seacas.patch
git apply builds/hip_stk.patch
cd -
