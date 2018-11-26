################################################################################
#
# Set up env on waterman for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU
fi

if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9
  else
    echo
    echo "***"
    echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option"
    echo "*** for the compiler GNU.  Replace '$ATDM_CONFIG_KOKKOS_ARCH' in the"
    echo "*** job name with 'Power9'"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] ; then
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9,Volta70
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9,Volta70
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Volta70" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9,Volta70
  else
    echo
    echo "***"
    echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option"
    echo "*** for the CUDA compiler.  Replace '$ATDM_CONFIG_KOKKOS_ARCH' in the"
    echo "*** job name with one of the following options:"
    echo "***"
    echo "***   'Volta70' Power9 with Volta V-100 GPU (Default)"
    echo "***"
    return
  fi
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return
fi

echo "Using waterman compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=64
# NOTE: Above settings are used for building on a Dual-Socket POWER9 with 8
# cores per socket.

module purge

module load git/2.10.1
# NOTE: Must load a git module since /usr/bin/git does not exist on the
# compute nodes.

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=20
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=40
fi

if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
    module load devpack/20180517/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88
    module swap openblas/0.2.20/gcc/7.2.0 netlib/3.8.0/gcc/7.2.0
    export OMPI_CXX=`which g++`
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] ; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-9.2  # The default CUDA version currently
  fi
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2" ]] ; then
    module load devpack/20180517/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88
    module swap openblas/0.2.20/gcc/7.2.0 netlib/3.8.0/gcc/7.2.0
  else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not a supported version of CUDA on this system!"
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
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"
  export ATDM_CONFIG_USE_CUDA=ON
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  # Avoids timeouts due to not running on separate GPUs (see #2446)
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

# Use manually installed cmake and ninja to try to avoid module loading
# problems (see TRIL-208).  NOTE: These were build for Power8 on 'white' and
# 'ride' but they also seem to work just fine on Power9 'waterman' as well :-)
export PATH=/ascldap/users/rabartl/install/white-ride/cmake-3.11.2/bin:/ascldap/users/rabartl/install/white-ride/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAGS="-map-by;socket:PE=4"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
