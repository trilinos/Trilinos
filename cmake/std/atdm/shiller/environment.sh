################################################################################
#
# Set up env on shiller/hansen for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU
fi

echo "Using hansen/shiller compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=32

module purge

if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
  export ATDM_CONFIG_KOKKOS_ARCH=HSW
  module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
  export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
  if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
    export OMP_NUM_THREADS=2
    # NOTE: For some reason, tests run with the GNU OPenMP builds don't seem
    # to have too much of a problem with running on all 16 cores with two
    # threads per core (with hyperthreading turned on).
  else
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
    # NOTE: For some reason, the tests run with the GNU serial builds seem to
    # be running on top of each other really badly so we need to reduce the
    # parallel level to avoid timeouts.  See #2976.
  fi
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
  module load devpack/openmpi/2.1.1/intel/17.4.196/cuda/none
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export ATDM_CONFIG_LAPACK_LIB="-mkl"
  export ATDM_CONFIG_BLAS_LIB="-mkl"
  if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
    export OMP_NUM_THREADS=2
  else
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  fi
  # NOTE: For some reason the tests run with the Intel builds (serial or
  # OpenMP) don't seem to have the problem of running on top of each other and
  # causing the tests to take longer.  Therefore, we can use all 16 cores (and
  # two threads per core with OpenMP and hyper threading).
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-8.0  # The default CUDA version currently
  fi
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-8.0" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.0" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/9.0.176
  else
      echo "***"
      echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not a supported version of CUDA on this system!"
      echo "***"
      return
  fi
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

# Use manually installed cmake and ninja (see TRIL-209)
export PATH=/home/rabartl/install/hansen-shiller/cmake-3.11.2/bin:/home/rabartl/install/hansen-shiller/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=16;--oversubscribe"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
