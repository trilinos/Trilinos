################################################################################
#
# Set up env on a SEMS NFS mounted RHEL6 for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

echo "Using SEMS RHEL6 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=OFF
export ATDM_CONFIG_BUILD_COUNT=32

module purge
module load sems-env
module load sems-cmake/3.5.2
module load sems-git/2.10.1

#module load ninja/1.7.2

if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  export OMP_NUM_THREADS=2
  # NOTE: With hyper-threading enabled, the 16 cores can run two threads each
  # or 32 threads total.
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  # NOTE: When running in serial, the second hyperthread can't see to run a
  # sperate MPI process and if you try to run with -j32 with ctest, you get a
  # bunch of failures that say "libgomp: Thread creation failed: Resource
  # temporarily unavailable".  Therefore, we need to run with -j16, not -j32.
fi

if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
    module load sems-gcc/6.1.0
    export OMPI_CXX=`which g++`
    export LAPACK_ROOT=/usr/lib64/atlas
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT};-llapack"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas"
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load sems-gcc/4.9.3
    module load sems-intel/17.0.1
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SEMS_INTEL_ROOT/mkl/lib/intel64/
    export OMPI_CXX=`which icpc`
    export OMPI_CC=`which icc`
    export OMPI_FC=`which ifort`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
else
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

module load sems-openmpi/1.10.1
module load sems-netcdf/4.4.1/exo_parallel
module load sems-hdf5/1.8.12/parallel
module load sems-zlib/1.2.8/base
module load sems-boost/1.59.0/base
module load sems-superlu/4.3/base

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=-lhwloc

export BOOST_ROOT=${SEMS_BOOST_ROOT}
export HDF5_ROOT=${SEMS_HDF5_ROOT}
export NETCDF_ROOT=${SEMS_NETCDF_ROOT}

export ATDM_CONFIG_HDF5_LIBS="-L${SEMS_HDF5_ROOT}/lib;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${SEMS_BOOST_ROOT}/lib;-L${SEMS_NETCDF_ROOT}/lib;-L${SEMS_NETCDF_ROOT}/lib;-L${SEMS_PNETCDF_ROOT}/lib;-L${SEMS_HDF5_ROOT}/lib;${SEMS_BOOST_ROOT}/lib/libboost_program_options.a;${SEMS_BOOST_ROOT}/lib/libboost_system.a;${SEMS_NETCDF_ROOT}/lib/libnetcdf.a;${SEMS_NETCDF_ROOT}/lib/libpnetcdf.a;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
