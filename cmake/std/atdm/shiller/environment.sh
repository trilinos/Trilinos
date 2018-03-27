################################################################################
#
# Set up env on shiller/hansen for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

echo "Using hansen/shiller compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=32

module purge

module load ninja/1.7.2

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
    export ATDM_CONFIG_KOKKOS_ARCH=HSW
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=`which g++`
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load devpack/openmpi/2.1.1/intel/17.4.196/cuda/none
    export OMPI_CXX=`which icpc`
    export OMPI_CC=`which icc`
    export OMPI_FC=`which ifort`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
elif [ "$ATDM_CONFIG_COMPILER" == "CUDA" ]; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=$ATDM_CONFIG_TRILNOS_DIR/packages/kokkos/bin/nvcc_wrapper 
    if [ ! -x "$OMPI_CXX" ]; then
        echo "No nvcc_wrapper found"
        return
    fi
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_USE_CUDA=ON
    export CUDA_LAUNCH_BLOCKING=1
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
    # Avoids timeouts due to not running on seprate GPUs (see #2446)
else
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;-lhdf5_hl;-lhdf5;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=16;--oversubscribe"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
