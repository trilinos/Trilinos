################################################################################
#
# Set up env on shiller/hansen for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU
fi

if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
  else
    export ATDM_CONFIG_KOKKOS_ARCH=HSW
  fi
fi

echo "Using hansen/shiller compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=32

module purge

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

if [[ "$ATDM_CONFIG_COMPILER" == "GNU" && "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW" ]] ; then
  module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL" && "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW" ]] ; then
  module load devpack/openmpi/2.1.1/intel/17.4.196/cuda/none
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export ATDM_CONFIG_LAPACK_LIBS="-mkl"
  export ATDM_CONFIG_BLAS_LIBS="-mkl"
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* && "$ATDM_CONFIG_KOKKOS_ARCH" == "Kepler37" ]] ; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-8.0  # The default CUDA version currently
  fi
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-8.0" ]] ; then
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.0" ]] ; then
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/9.0.176
  else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not a supported version"
    echo "*** of CUDA on this system!  Supported versions include 'cuda-8.0'"
    echo "*** and 'cuda-9.0'"
    echo "***"
    return
  fi
  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER} 
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_USE_CUDA=ON
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  # Avoids timeouts due to not running on seprate GPUs (see #2446)
else
  echo
  echo "***"
  echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH with compiler=$ATDM_CONFIG_COMPILER "
  echo "*** is not a supported combination on this system. Supported compilers are: "
  echo "*** gnu, intel, cuda-8.0, and cuda-9.0"
  echo "*** Please be aware that you do not need to specify KOKKOS_ARCH on this system because the "
  echo "*** the defaults are the only reasonable choices"
  echo "***"
  echo "*** FYI:"
  echo "*** 'HSW' (Haswell) will be used for all non-CUDA compiler builds and "
  echo "*** 'Kepler37' (Kepler K-80 GPU) will be used for all CUDA builds"
  echo "***"
  return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;-lhdf5_hl;-lhdf5;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"

# Use manually installed cmake and ninja (see TRIL-209)
export PATH=/home/rabartl/install/hansen-shiller/cmake-3.11.2/bin:/home/rabartl/install/hansen-shiller/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAGS="-map-by;socket:PE=16;--oversubscribe"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
