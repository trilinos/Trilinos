################################################################################
#
# Set up env on a SEMS NFS mounted RHEL7 for ATMD builds of Trilinos
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
    export ATDM_CONFIG_COMPILER=CLANG-7.0.1
  elif [[ "$ATDM_CONFIG_COMPILER" != "CLANG-7.0.1" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: CLANG COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only CLANG compilers supported on this system are:"
    echo "***   clang (defaults to clang-7.0.1)"
    echo "***   clang-7.0.1"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "GNU" ]] ; then
    export ATDM_CONFIG_COMPILER=GNU-7.2.0
  elif [[ "$ATDM_CONFIG_COMPILER" != "GNU-7.2.0" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: GNU COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only GNU compilers supported on this system are:"
    echo "***   gnu (defaults to gnu-7.2.0)"
    echo "***   gnu-7.2.0 (default)"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "INTEL" ]] ; then
    export ATDM_CONFIG_COMPILER=INTEL-18.0.5
  elif [[ "$ATDM_CONFIG_COMPILER" != "INTEL-17.0.1" ]] \
    && [[ "$ATDM_CONFIG_COMPILER" != "INTEL-18.0.5" ]]; then
    echo
    echo "***"
    echo "*** ERROR: INTEL COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only INTEL compilers supported on this system are:"
    echo "***   intel (defaults to intel-18.0.5)"
    echo "***   intel-17.0.1"
    echo "***   intel-18.0.5"
    echo "***"
    return
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] ; then
    export ATDM_CONFIG_COMPILER=CUDA-9.2
  elif [[ "$ATDM_CONFIG_COMPILER" != "CUDA-9.2" ]] ; then
    echo
    echo "***"
    echo "*** ERROR: CUDA COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
    echo "*** Only CUDA compilers supported on this system are:"
    echo "***   cuda (defaults to cuda-9.2)"
    echo "***   cuda-9.2"
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
  unset ATDM_CONFIG_KOKKOS_ARCH
fi

echo "Using SEMS RHEL7 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

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

if [[ "${SEMS_MODULEFILES_ROOT}" == "" ]] ; then
  if [[ -d /projects/sems/modulefiles ]] ; then
    echo "NOTE: SEMS modules not defined, loading their definition!"
    module use /projects/sems/modulefiles/projects
    export SEMS_MODULEFILES_ROOT=/projects/sems/modulefiles
  else
    echo "ERROR: The SEMS modules are not defined and default location does not exist!"
    return
  fi
fi

module purge
module load sems-env
module load sems-git/2.10.1

module load sems-cmake/3.12.2
module load sems-ninja_fortran/1.8.2

if [[ "$ATDM_CONFIG_NODE_TYPE" == "CUDA" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=4
  # We just need to be super conservative by default when using a GPU.  If
  # users want to use more MPI processes, then can override this with
  # ATDM_CONFIG_CTEST_PARALLEL_LEVEL_OVERRIDE.
elif [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
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

if [[ "$ATDM_CONFIG_COMPILER" == "CLANG-7.0.1" ]] ; then
  module load sems-clang/7.0.1
  export OMPI_CXX=`which clang++`
  export OMPI_CC=`which clang`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="/usr/lib64/liblapack.so.3"
  export ATDM_CONFIG_BLAS_LIBS="/usr/lib64/libblas.so.3"
elif [[ "$ATDM_CONFIG_COMPILER" == "CLANG-3.9.0" ]] ; then
  module load sems-clang/3.9.0
  export OMPI_CXX=`which clang++`
  export OMPI_CC=`which clang`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="/usr/lib64/liblapack.so.3"
  export ATDM_CONFIG_BLAS_LIBS="/usr/lib64/libblas.so.3"
elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0" ]] ; then
  module load sems-gcc/7.2.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="/usr/lib64/liblapack.so.3"
  export ATDM_CONFIG_BLAS_LIBS="/usr/lib64/libblas.so.3"
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL-17.0.1" ]] ; then
  module load sems-intel/17.0.1
  module load atdm-env
  module load atdm-mkl/18.0.5
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export ATDM_CONFIG_LAPACK_LIBS="-mkl"
  export ATDM_CONFIG_BLAS_LIBS="-mkl"
  export LM_LICENSE_FILE=28518@cee-infra009.sandia.gov
  if [[ "${ATDM_CONFIG_LM_LICENSE_FILE_OVERRIDE}" != "" ]] ; then
    export LM_LICENSE_FILE=${ATDM_CONFIG_LM_LICENSE_FILE_OVERRIDE}
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL-18.0.5" ]] ; then
  module load sems-gcc/7.2.0
  module load sems-intel/18.0.5
  module load atdm-env
  module load atdm-mkl/18.0.5
  export ATDM_CONFIG_CXX_FLAGS="-D_GLIBCXX_USE_CXX11_ABI=0"
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export ATDM_CONFIG_LAPACK_LIBS="-mkl"
  export ATDM_CONFIG_BLAS_LIBS="-mkl"
  export LM_LICENSE_FILE=28518@cee-infra009.sandia.gov
  if [[ "${ATDM_CONFIG_LM_LICENSE_FILE_OVERRIDE}" != "" ]] ; then
    export LM_LICENSE_FILE=${ATDM_CONFIG_LM_LICENSE_FILE_OVERRIDE}
  fi
elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA-9.2" ]] ; then
  module load sems-gcc/7.2.0
  module load sems-cuda/9.2
  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER}
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export ATDM_CONFIG_LAPACK_LIBS="/usr/lib64/liblapack.so.3"
  export ATDM_CONFIG_BLAS_LIBS="/usr/lib64/libblas.so.3"
  # some Trilinos tests require this to run correctly
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  export KOKKOS_NUM_DEVICES=2
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

module load sems-openmpi/1.10.1
module load sems-netcdf/4.4.1/exo_parallel
module load sems-hdf5/1.8.12/parallel
module load sems-zlib/1.2.8/base
module load sems-boost/1.59.0/base
module unload sems-python/2.7.9 
module load sems-superlu/4.3/base

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] && \
  [[ "${ATDM_CONFIG_COMPLEX}" == "ON" ]] && \
  [[ "${ATDM_CONFIG_SHARED_LIBS}" == "ON" ]] ; then
  export ATDM_CONFIG_USE_NINJA=OFF
  export ATDM_CONFIG_CMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS=OFF
fi
# NOTE: The reason for the above logic is that 'nvcc' can't handle *.rsp
# response files that CMake switches to when the command-lines become too long
# and there is no way to turn off the usage of response files with the CMake
# Ninja generator as of CMake 3.14.0.  The problem is that when CUDA and
# complex are enabled and shared libs are used, 'nvcc' is used to create the
# kokkoskernels shared lib which has a ton of object files which triggers
# internal CMake logic to condense these down into a single kokkoskernels.rsp
# file that it passes to 'nvcc' to create the kokkoskernels library.  When the
# CMake-generated *.rsp file is passed to 'nvcc', it does not know how to
# process it and it gives the error "No input files specified".  The
# workaround is to switch to the CMake Makefile generator which allows and
# turning off the usage of response files for the list of object files (and
# 'nvcc' seems to be able to handle this long command-line in this case).
# Note that we don't yet need to switch to use the CMake Makefile generator
# for 'static' builds since 'ar' is used to create the kokkoskernels lib which
# does not seem to have a problem with the creation of this lib.  (See
# TRIL-255 and TRIL-264.)

if [[ "${ATDM_CONFIG_SHARED_LIBS}" == "ON" ]] ; then
  ATDM_CONFIG_TPL_LIB_EXT=so
else
  ATDM_CONFIG_TPL_LIB_EXT=a
fi

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=-lhwloc

export BOOST_ROOT=${SEMS_BOOST_ROOT}
export HDF5_ROOT=${SEMS_HDF5_ROOT}
export NETCDF_ROOT=${SEMS_NETCDF_ROOT}

export ATDM_CONFIG_HDF5_LIBS="-L${SEMS_HDF5_ROOT}/lib;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${SEMS_BOOST_ROOT}/lib;-L${SEMS_NETCDF_ROOT}/lib;-L${SEMS_NETCDF_ROOT}/lib;-L${SEMS_PNETCDF_ROOT}/lib;-L${SEMS_HDF5_ROOT}/lib;${SEMS_BOOST_ROOT}/lib/libboost_program_options.${ATDM_CONFIG_TPL_LIB_EXT};${SEMS_BOOST_ROOT}/lib/libboost_system.${ATDM_CONFIG_TPL_LIB_EXT};${SEMS_NETCDF_ROOT}/lib/libnetcdf.a;${SEMS_NETCDF_ROOT}/lib/libpnetcdf.a;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

# NOTE: SEMS does not provide the correct *.so files for NetCDF so we can't
# use them in a shared lib build :-(

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
