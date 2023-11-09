################################################################################
#
# Set up env on Blake for builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

#
# Deal with compiler versions
#

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL-2021.4.0
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "INTEL" ]] ; then
    export ATDM_CONFIG_COMPILER=INTEL-2021.4.0
  elif [[ "$ATDM_CONFIG_COMPILER" != "INTEL-2021.4.0" ]]; then
    echo
    echo "***"
    echo "*** ERROR: INTEL COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only INTEL compilers supported on this system are:"
    echo "***   intel (defaults to intel-2021.4.0)"
    echo "***   intel-2021.4.0"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "GNU" ]] ; then
    export ATDM_CONFIG_COMPILER=GNU-11.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" != "GNU-11.2.0" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: GNU COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only GNU compilers supported on this system are:"
    echo "***   gnu (defaults to gnu-11.2.0)"
    echo "***   gnu-11.2.0 (default)"
    echo "***"
    return
  fi
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

#
# Allow KOKKOS_ARCH which is needed for CUDA builds
#

if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=SKX
elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "SKX" ]] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=SKX
else
  echo
  echo "***"
  echo "*** Unrecognized or unsupported Kokkos arch on this system - unsetting ATDM_CONFIG_KOKKOS_ARCH and letting Kokkos handle it"
  echo "***"
  unset ATDM_CONFIG_KOKKOS_ARCH
fi

echo "Using blake compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

# TODO: Change to off if issues with ninja...
export ATDM_CONFIG_USE_NINJA=ON

# Get ATDM_CONFIG_NUM_CORES_ON_MACHINE for this machine
source $ATDM_SCRIPT_DIR/utils/get_num_cores_on_machine.sh

if [ "$ATDM_CONFIG_NUM_CORES_ON_MACHINE" -gt "20" ] ; then
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=20
  # NOTE: We get links crashing if we try to use to many processes.  ToDo: We
  # should limit the number of processes that ninja uses to link instead of
  # reducing the overrall parallel build level like this.
else
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=$ATDM_CONFIG_NUM_CORES_ON_MACHINE
fi

export ATDM_CONFIG_BUILD_COUNT=$ATDM_CONFIG_MAX_NUM_CORES_TO_USE
# NOTE: Use as many build processes and there are cores by default.


module purge
module load cmake/3.25.2 git

if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  export OMP_NUM_THREADS=2
  # NOTE: With hyper-threading enabled, you can run as many threads as there
  # are cores and with 2 OpenMP threads per MPI process, the means you can run
  # as many MPI processes as there are cores on the machine with 2 OpenMP
  # threads.  But we want to be conservative and instead run with half that
  # many to be safe and avoid time-outs.
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  # NOTE: NOTE: When running in serial, the second hyperthread can't seem to
  # run a sperate MPI process and if you try to run with as many they you get
  # a bunch of failures that say "libgomp: Thread creation failed: Resource
  # temporarily unavailable".  So we can only run with as many MPI processes
  # as there are cores on the machine.  But we want to be conservative and
  # instead run with half that many to be safe and avoid time-outs.
  export OMP_PROC_BIND=FALSE
  export OMP_NUM_THREADS=1
fi

if [[ "$ATDM_CONFIG_COMPILER" == "INTEL-2021.4.0" ]] ; then
  module load intel/oneAPI/hpc-toolkit/2021.4.0 intel/oneAPI/base-toolkit/2021.4.0 
  module load openmpi/4.0.5/intel-oneapi/2021.4.0 
  module swap gcc/7.2.0 gcc/10.2.0
  module load openblas/0.3.21/gcc/10.2.0
  module load boost/1.75.0/intel-oneapi/2021.2.0
  module load hdf5/1.10.7/openmpi/4.0.5/intel-oneapi/2021.2.0 netcdf-c/4.7.4/openmpi/4.0.5/intel-oneapi/2021.2.0 zlib/1.2.11
  module load python/3.7.3
  #module load python/2.7.13

  export ATDM_CONFIG_CXX_FLAGS="-g"
  export ATDM_CONFIG_C_FLAGS=" -g"
  export ATDM_CONFIG_Fortran_FLAGS=" -g"
  #export ATDM_CONFIG_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-lopenblas"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lopenblas"
  export ATDM_CONFIG_MPI_POST_FLAGS="-bind-to;socket;-map-by;socket"
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-11.2.0" ]] ; then
  module load gcc/11.2.0 
  module load openmpi/3.1.6/gcc/11.2.0 
  module load netlib-lapack/3.10.1/gcc/11.2.0 
  module load boost/1.81.0/gcc/11.2.0
  module load hdf5/1.14.0/gcc/11.2.0/openmpi/3.1.6 zlib/1.2.13/gcc/11.2.0 netcdf-c/4.9.0/gcc/11.2.0/openmpi/3.1.6 netcdf-cxx/4.2/gcc/11.2.0/openmpi/3.1.6

  module load superlu/5.3.0/gcc/11.2.0 metis/5.1.0/gcc/11.2.0
  module load python/3.7.3
  #module load python/2.7.13

  export ATDM_CONFIG_CXX_FLAGS="-g"
  export ATDM_CONFIG_C_FLAGS="-g"
  export ATDM_CONFIG_Fortran_FLAGS="-g"
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  #export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp;-lm"
  #export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"
  export ATDM_CONFIG_BLAS_LIBS="${BLAS_ROOT}/lib64/libblas.so"
  export ATDM_CONFIG_LAPACK_LIBS="${LAPACK_ROOT}/lib64/liblapack.so"
  export ATDM_CONFIG_MPI_POST_FLAGS="-bind-to;socket;-map-by;socket"
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

if [[ "${ATDM_CONFIG_SHARED_LIBS}" == "ON" ]] ; then
  ATDM_CONFIG_TPL_LIB_EXT=so
else
  ATDM_CONFIG_TPL_LIB_EXT=a
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;${ZLIB_ROOT}/lib/libz.a;-ldl"

if [[ "${NETCDF_C_ROOT}" ]] ; then
  export NETCDF_ROOT=${NETCDF_C_ROOT}
fi

if [[ "${PARALLEL_NETCDF_ROOT}" == "" ]] ; then
  export PNETCDF_ROOT=${NETCDF_ROOT}
  export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib64;${NETCDF_ROOT}/lib64/libnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"
else
  export NETCDF_ROOT=${NETCDF_ROOT}
  export PNETCDF_ROOT=${PARALLEL_NETCDF_ROOT}
  export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib64;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"
fi

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif77`

export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
