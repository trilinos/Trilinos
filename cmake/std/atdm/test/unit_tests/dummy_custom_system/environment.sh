# Custom configuration used for unit testing

echo "Using dummy_custom_system for compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=5
export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=3
unset ATDM_CONFIG_KOKKOS_ARCH # Or set to the arch you know you have!

# MPI Stuff
export MPICC=dummy_mpicc
export MPICXX=dummy_mpicxx
export MPIF90=dummy_mpif90
export OMPI_CXX=dummy_g++
export OMPI_CC=dummy_gcc
export OMPI_FC=dummy_gfortarn
export ATDM_CONFIG_MPI_PRE_FLAGS="-pre;--flag=2;-args=1"

# OpenMP stuff
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
unset OMP_PLACES

# TPL stuff
export ATDM_CONFIG_USE_HWLOC=OFF
export LAPACK_ROOT=/usr/lib64
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-llapack"
export ATDM_CONFIG_BLAS_LIBS="-L${LAPACK_ROOT}/lib;-lblas"
export ZLIB_ROOT=/lib/to/zlib
export BOOST_ROOT=/lib/to/boost
export HDF5_ROOT=/lib/to/hdf5
export NETCDF_ROOT=/lib/to/netcdf

# Done
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
