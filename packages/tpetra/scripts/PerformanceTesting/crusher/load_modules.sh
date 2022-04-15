#!/bin/bash
module load craype-accel-amd-gfx90a
module load cray-libsci/21.08.1.2
module load rocm/4.5.2
module load cray-mpich
module load cray-hdf5/1.12.0.7 cray-netcdf/4.7.4.7
module load parallel-netcdf/1.12.2
module load cmake/3.22.1
module load boost/1.78.0
module load zlib/1.2.11
module load parmetis/4.0.3
module load metis/5.1.0
module load yaml-cpp/0.7.0
module load zlib/1.2.11
module load superlu/5.2.2
module load ninja

export MPICH_GPU_SUPPORT_ENABLED=1

#If you compile with Cray's cc and CC, the following
#must be set before compiling so the executable picks up GTL.
PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
