################################################################################
#
# Set up env on a SEMS NFS mounted RHEL6 for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

#
# Deal with compiler versions
#

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU-7.2.0
elif [[ "$ATDM_CONFIG_COMPILER" == "CLANG"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CLANG" ]] ; then
    export ATDM_CONFIG_COMPILER=CLANG-3.9.0
  elif [[ "$ATDM_CONFIG_COMPILER" != "CLANG-3.9.0" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: CLANG COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only CLANG compilers supported on this system are:"
    echo "***   clang (defaults to clang-3.9.0)"
    echo "***   clang-3.9.0"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "GNU" ]] ; then
    export ATDM_CONFIG_COMPILER=GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" != "GNU-6.1.0" ]] &&
    [[ "$ATDM_CONFIG_COMPILER" != "GNU-7.2.0" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: GNU COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only GNU compilers supported on this system are:"
    echo "***   gnu (defaults to gnu-7.2.0)"
    echo "***   gnu-6.1.0"
    echo "***   gnu-7.2.0 (default)"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "INTEL" ]] ; then
    export ATDM_CONFIG_COMPILER=INTEL-17.0.1
  elif [[ "$ATDM_CONFIG_COMPILER" != "INTEL-17.0.1" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: INTEL COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only INTEL compilers supported on this system are:"
    echo "***   intel (defaults to intel-17.0.1)"
    echo "***   intel-17.0.1"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then
  echo
  echo "***"
  echo "*** ERROR: CUDA COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "*** The 'sems-rhel6' env does not support CUDA.  NOTE: You might want"
  echo "*** the 'sems-rhel6' env that does support CUDA."
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
  unset ATDM_CONFIG_KOKKOS_ARCH
fi

echo "Using SEMS RHEL6 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

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
module load sems-env
module load sems-git/2.10.1

module load sems-cmake/3.12.2
module load sems-ninja_fortran/1.8.2

if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  export OMP_NUM_THREADS=2
  export OMP_PROC_BIND=false
  unset OMP_PLACES
  # NOTE: With hyper-threading enabled, you can run as many threads as there
  # are cores and with 2 OpenMP threads per MPI process, the means you can run
  # as many MPI processes as there are cores on the machine with 2 OpenMP
  # threads.  But we want to be conservative and instead run with half that
  # many to be safe and avoid time-outs.
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  export OMP_NUM_THREADS=1
  export OMP_PROC_BIND=false
  unset OMP_PLACES
  # NOTE: NOTE: When running in serial, the second hyperthread can't seem to
  # run a sperate MPI process and if you try to run with as many they you get
  # a bunch of failures that say "libgomp: Thread creation failed: Resource
  # temporarily unavailable".  So we can only run with as many MPI processes
  # as there are cores on the machine.  But we want to be conservative and
  # instead run with half that many to be safe and avoid time-outs.
fi

if [[ "$ATDM_CONFIG_COMPILER" == "CLANG-3.9.0" ]] ; then
  module load sems-clang/3.9.0
  export OMPI_CXX=`which clang++`
  export OMPI_CC=`which clang`
  export OMPI_FC=`which gfortran`
  export LAPACK_ROOT=/usr/lib64/atlas
  export BLAS_ROOT=/usr/lib64/atlas
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-llapack"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas"
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-6.1.0" ]] ; then
  module load sems-gcc/6.1.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export LAPACK_ROOT=/usr/lib64/atlas
  export BLAS_ROOT=/usr/lib64/atlas
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-llapack"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas"
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0" ]] ; then
  module load sems-gcc/7.2.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export LAPACK_ROOT=/usr/lib64/atlas
  export BLAS_ROOT=/usr/lib64/atlas
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-llapack"
  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas"
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL-17.0.1" ]] ; then
  module load sems-intel/17.0.1
  module swap sems-gcc/4.8.4 sems-gcc/4.9.3
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SEMS_INTEL_ROOT/mkl/lib/intel64/
  export ATDM_CONFIG_LAPACK_LIBS="-mkl"
  export ATDM_CONFIG_BLAS_LIBS="-mkl"
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

module load sems-openmpi/1.10.1
module load sems-netcdf/4.7.3/parallel
module load sems-hdf5/1.10.6/parallel
module load sems-zlib/1.2.8/base
module unload sems-python/2.7.9
module load sems-superlu/4.3/base

if [[ "${ATDM_CONFIG_SHARED_LIBS}" == "ON" ]] ; then
  ATDM_CONFIG_TPL_LIB_EXT=so
else
  ATDM_CONFIG_TPL_LIB_EXT=a
fi

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=-lhwloc

export ZLIB_ROOT=${SEMS_ZLIB_ROOT}
export HDF5_ROOT=${SEMS_HDF5_ROOT}
export NETCDF_ROOT=${SEMS_NETCDF_ROOT}
if [[ "${SEMS_PNETCDF_ROOT}" == "" ]] ; then
  export PNETCDF_ROOT=${SEMS_NETCDF_ROOT}
else
  export PNETCDF_ROOT=${SEMS_PNETCDF_ROOT}
fi

export ATDM_CONFIG_HDF5_LIBS="${HDF5_ROOT}/lib/libhdf5_hl.${ATDM_CONFIG_TPL_LIB_EXT};${HDF5_ROOT}/lib/libhdf5.${ATDM_CONFIG_TPL_LIB_EXT};${ZLIB_ROOT}/lib/libz.${ATDM_CONFIG_TPL_LIB_EXT};-ldl"

export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;${NETCDF_ROOT}/lib/libnetcdf.${ATDM_CONFIG_TPL_LIB_EXT};${PNETCDF_ROOT}/lib/libpnetcdf.${ATDM_CONFIG_TPL_LIB_EXT};${ATDM_CONFIG_HDF5_LIBS};-lcurl"

# NOTE: SEMS does not provide a *.a files for PNetCDF so we can't use them in
# a shared lib build :-(

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

#
# Set up default install-related stuff
#

export ATDM_CONFIG_WORKSPACE_BASE_DEFAULT=/home/atdm-devops-admin/jenkins
export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=/home/atdm-devops-admin/trilinos_installs
export ATDM_CONFIG_INSTALL_PBP_RUNNER_DEFAULT=/home/atdm-devops-admin/tools/run-as-atdm-devops-admin

#
# Done!
#

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
