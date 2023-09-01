#!/bin/bash
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load craype-x86-trento
#module load cray-libsci  #loaded by PrgEnv-cray
module load rocm/5.2.0
module load cmake/3.22.1
module load ninja

module load cray-hdf5-parallel/1.12.1.1 cray-netcdf-hdf5parallel
module load parallel-netcdf/1.12.2
module load boost/1.79.0
module load zlib/1.2.11
#module load parmetis/4.0.3
#module load metis/5.1.0
#module load yaml-cpp/0.7.0
module load zlib/1.2.11
module load superlu/5.3.0

export MPICH_GPU_SUPPORT_ENABLED=1

#If you compile with Cray's cc and CC, the following
#must be set before compiling so the executable picks up GTL.
#PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
#PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
