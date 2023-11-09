#!/bin/bash

# Error out the job if any command here fails
set -e

source $WORKSPACE/PerfScripts/loadEnv.sh

cd $WORKSPACE/build

export OMPI_CXX=${TRILINOS_SRC}/packages/kokkos/bin/nvcc_wrapper

rm -rf CMakeFiles
rm -f CMakeCache.txt

cmake \
-DCMAKE_C_COMPILER=mpicc \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_BUILD_TYPE=Release \
-DTrilinos_TEST_CATEGORIES="PERFORMANCE" \
-DTrilinos_ENABLE_TESTS=ON \
-DTrilinos_ENABLE_EXAMPLES=ON \
-DTrilinos_ENABLE_Tpetra=ON \
-DTrilinos_ENABLE_Ifpack2=ON \
-DTrilinos_ENABLE_Belos=ON \
-DTrilinos_ENABLE_MueLu=ON \
-DTrilinos_ENABLE_PanzerMiniEM=ON \
-DTrilinos_ENABLE_Percept=OFF \
-DTrilinos_ENABLE_Stokhos=OFF \
-DTpetra_INST_INT_INT=ON \
-DXpetra_ENABLE_Epetra=ON \
-DMueLu_ENABLE_Epetra=ON \
-DTPL_ENABLE_MPI=ON \
-DTPL_ENABLE_HDF5=ON \
-DTPL_ENABLE_Netcdf=ON \
-DTPL_ENABLE_Matio=OFF \
-DTPL_ENABLE_X11=OFF \
-DTPL_ENABLE_BLAS=ON \
-DTPL_ENABLE_LAPACK=ON \
-DTPL_ENABLE_CUBLAS=ON \
-DTPL_ENABLE_CUSPARSE=ON \
-DTPL_LAPACK_LIBRARIES="-L/usr/tcetmp/packages/lapack/lapack-3.8.0-gcc-4.9.3/lib;-llapack;-lgfortran;-lgomp" \
-DTPL_BLAS_LIBRARIES="-L/lib;-lblas;-lgfortran;-lgomp;-lm" \
-DTPL_Netcdf_LIBRARIES="$NETCDF_ROOT/lib64/libnetcdf.a;$HDF5_ROOT/lib/libhdf5_hl.a;$HDF5_ROOT/lib/libhdf5.a;$PNETCDF_ROOT/lib/libpnetcdf.a" \
-DBoost_INCLUDE_DIRS=$BOOST_ROOT/include \
-DBoost_LIBRARY_DIRS=$BOOST_ROOT/lib \
-DNetcdf_INCLUDE_DIRS=$NETCDF_ROOT/include \
-DNetcdf_LIBRARY_DIRS=$NETCDF_ROOT/lib \
-DHDF5_INCLUDE_DIRS=$HDF5_ROOT/include \
-DHDF5_LIBRARY_DIRS=$HDF5_ROOT/lib \
-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
-DKokkos_ENABLE_CUDA=ON \
-DKokkos_ENABLE_CUDA_LAMBDA=ON \
-DKokkos_ENABLE_CUDA_CONSTEXPR=ON \
-DTpetra_INST_CUDA=ON \
-DTpetra_INST_SERIAL=ON \
-DKokkos_ARCH_POWER9=ON \
-DKokkos_ARCH_VOLTA70=ON \
-DTPL_ENABLE_CUDA:BOOL=ON \
-DTrilinos_AUTOGENERATE_TEST_RESOURCE_FILE=OFF \
-DTrilinos_ENABLE_PyTrilinos=OFF \
-DMPI_EXEC=$TRILINOS_SRC/cmake/std/atdm/ats2/trilinos_jsrun \
-DMPI_EXEC_NUMPROCS_FLAG="-p" \
-DMPI_EXEC_POST_NUMPROCS_FLAGS="--rs_per_socket;4" \
-DMPI_EXEC_MAX_NUMPROCS=16 \
$TRILINOS_SRC

make -j36

