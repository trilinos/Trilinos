################################################################################
#
# Set up env on ride/white for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

#
# Deal with compiler versions
#

if [[ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ]] ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "GNU" ]] ; then
    export ATDM_CONFIG_COMPILER=GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0" ]] ; then
    export ATDM_CONFIG_COMPILER=GNU-7.2.0
  else
    echo
    echo "***"
    echo "*** ERROR: GNU COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only GNU compilers supported on this system are:"
    echo "***   gnu (defaults to gnu-7.2.0)"
    echo "***   gnu-7.2.0 (default)"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-9.2_GNU-7.2.0  # The default CUDA version currently
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-9.2_GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2_GNU-7.2.0" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-9.2_GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1_GNU-7.2.0" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-10.1_GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1_GNU-7.2.0" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-10.1_GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-10.1_GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-10.1_GNU-7.2.0
  else
    echo
    echo "***"
    echo "*** ERROR: CUDA COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only CUDA compilers supported on this system are:"
    echo "***   cuda (defaults to cuda-9.2-gnu-7.2.0)"
    echo "***   cuda-9.2 (defaults to cuda-9.2-gnu-7.2.0)"
    echo "***   cuda-9.2-gnu-7.2.0 (default)"
    echo "***   cuda-10.1 (defaults to cuda-10.1-gnu-7.2.0)"
    echo "***   cuda-10.1-gnu-7.2.0"
    echo "***"
    return
  fi
fi


#
# Deal with KOKKOS_ARCH
#

if [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8
    export ATDM_CONFIG_QUEUE=rhel7F
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Power8" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8
    export ATDM_CONFIG_QUEUE=rhel7F
  else
    echo
    echo "***"
    echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option"
    echo "*** for the compiler $ATDM_CONFIG_COMPILER."
    echo "**  Replace '$ATDM_CONFIG_KOKKOS_ARCH' in '${ATDM_CONFIG_BUILD_NAME}'"
    echo "*** with 'Power8'!"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] ; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1_GNU-7.2.0" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,Pascal60
    export ATDM_CONFIG_QUEUE=dev  # ToDo: Update once ready for production!
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,Kepler37
    export ATDM_CONFIG_QUEUE=rhel7F
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Power8" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,Kepler37
    export ATDM_CONFIG_QUEUE=rhel7F
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Kepler37" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,$ATDM_CONFIG_KOKKOS_ARCH
    export ATDM_CONFIG_QUEUE=rhel7F
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Kepler35" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,$ATDM_CONFIG_KOKKOS_ARCH
    export ATDM_CONFIG_QUEUE=rhel7T
  elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "Pascal60" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power8,$ATDM_CONFIG_KOKKOS_ARCH
    export ATDM_CONFIG_QUEUE=rhel7G
  else
    echo
    echo "***"
    echo "*** ERROR: KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not a valid option"
    echo "*** for the compiler $ATDM_CONFIG_COMPILER."
    echo "**  Replace '$ATDM_CONFIG_KOKKOS_ARCH' in '${ATDM_CONFIG_BUILD_NAME}'"
    echo "*** with one of the following options:"
    echo "***"
    echo "***   'Kepler35' Power8 with Kepler K-40 GPU"
    echo "***   'Kepler37' Power8 with Kepler K-80 GPU (Default)"
    echo "***   'Pascal60' Power8 with Pascal P-100 GPU"
    echo "***"
    return
  fi
else
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return
fi

echo "Using white/ride compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

#
# Finish clearing the current env
#

module purge

#
# Deal with build, test and OpenMP runtime options
#

export ATDM_CONFIG_USE_NINJA=ON

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] && \
  [[ "${ATDM_CONFIG_CUDA_RDC}" == "ON" ]] ; then
  # When CUDA+RDC is enabled, it uses more RAM building with more build
  # processes crashes nodes on 'white' and 'ride (see #4502)
  if [[ "${ATDM_CONFIG_PT_PACKAGES}" == "ON" ]] ; then
    # Building with all PT packages crashes nodes on 'white' and 'ride' for
    # the same parallel levels that work for just the ATDM packages below.
    # See ATDV-155.
    export ATDM_CONFIG_BUILD_COUNT=20
    export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=10
  else
    # We seem to be able to use more build processes when just ATDM packages
    # are enabled.
    export ATDM_CONFIG_BUILD_COUNT=32
    export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=16
  fi
else
  export ATDM_CONFIG_BUILD_COUNT=64
fi

# NOTE: Above settings are used for running on a single rhel7F (Firestone,
# Dual-Socket POWER8, 8 cores per socket, K80 GPUs) node.

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=32
fi

if [ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0" ] ; then

  # Load the modules and set up env
  module load devpack/20180521/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88
  module swap openblas/0.2.20/gcc/7.2.0 netlib/3.8.0/gcc/7.2.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"

elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2_GNU-7.2.0" ]] ; then

  module load devpack/20180521/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88
  module swap openblas/0.2.20/gcc/7.2.0 netlib/3.8.0/gcc/7.2.0
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
  export KOKKOS_NUM_DEVICES=2
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=4
  # ABove avoids timeouts due to not running on separate GPUs (see #2446)

elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1_GNU-7.2.0" ]] ; then

  module load devpack/20190404/openmpi/4.0.1/gcc/7.2.0/cuda/10.1.105
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
  export KOKKOS_NUM_DEVICES=2
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=4
  # ABove avoids timeouts due to not running on separate GPUs (see #2446)

fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

# Use manually installed cmake and ninja to try to avoid module loading
# problems (see TRIL-208)
export PATH=/ascldap/users/rabartl/install/white-ride/cmake-3.11.2/bin:/ascldap/users/rabartl/install/white-ride/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAGS="-map-by;socket:PE=4"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
