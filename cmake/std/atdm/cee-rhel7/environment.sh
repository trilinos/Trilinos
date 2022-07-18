################################################################################
#
# Set up env on a CEE RHEL7 system for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

# NOTE: The custom_builds.sh script in this directory sets the exact compilers
# used so no need for dealing with different varients of compilers here.

if [[ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ]] ; then
  # Abort, no compiler was selected!
  return
fi

if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] ; then
  unset ATDM_CONFIG_KOKKOS_ARCH
  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1.243_GNU-7.2.0_OPENMPI-4.0.3" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=VOLTA70
  fi
else
  echo
  echo "***"
  echo "*** ERROR: Specifying KOKKOS_ARCH is not supported on CEE RHEL7 builds"
  echo "*** remove '$ATDM_CONFIG_KOKKOS_ARCH' from JOB_NAME=$JOB_NAME"
  echo "***"
  return
fi

echo "Using CEE RHEL7 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

# Get ATDM_CONFIG_NUM_CORES_ON_MACHINE for this machine
source $ATDM_SCRIPT_DIR/utils/get_num_cores_on_machine.sh

if [ "$ATDM_CONFIG_NUM_CORES_ON_MACHINE" -gt "16" ] ; then
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=16
  # NOTE: We get links crashing if we try to use to many processes.  ToDo: We
  # should limit the number of processes that ninja uses to link instead of
  # reducing the overrall parallel build level like this.
else
  export ATDM_CONFIG_MAX_NUM_CORES_TO_USE=$ATDM_CONFIG_NUM_CORES_ON_MACHINE
fi

export ATDM_CONFIG_BUILD_COUNT=$ATDM_CONFIG_MAX_NUM_CORES_TO_USE
# NOTE: Use as many build processes and there are cores by default.

if [[ "${ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE}" != "1" ]] ; then
  module purge
else
  echo "NOTE: ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE=1 is set so using pre-loaded sparc-dev module!"
fi

# Warning options requested by Gemma team (which should hopefully also take
# care of warnings required by the other ATDM APPs as well).  See #3178 and
# #4221
ATDM_CONFIG_GNU_CXX_WARNINGS="-Wall -Wextra"
ATDM_CONFIG_INTEL_CXX_WARNINGS="-Wall -Warray-bounds -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Wuninitialized -Wmaybe-uninitialized -Wmain -Wnarrowing -Wnonnull -Wparentheses -Wpointer-sign -Wreorder -Wreturn-type -Wsign-compare -Wsequence-point -Wtrigraphs -Wunused-function -Wunused-but-set-variable -Wunused-variable -Wwrite-strings"

# For now, turn on warnings by default:
if [[ "${ATDM_CONFIG_ENABLE_STRONG_WARNINGS}" == "" ]] ; then
  export ATDM_CONFIG_ENABLE_STRONG_WARNINGS=1
fi

if  [[ "$ATDM_CONFIG_COMPILER" == "CLANG-9.0.1_OPENMPI-4.0.3" ]]; then
  atdm_config_load_sparc_dev_module sparc-dev/clang-9.0.1_openmpi-4.0.3
  export OMPI_CXX=`which clang++`
  export OMPI_CC=`which clang`
  export OMPI_FC=`which gfortran`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  if [[ "$ATDM_CONFIG_ENABLE_STRONG_WARNINGS" == "1" ]]; then
    export ATDM_CONFIG_CXX_FLAGS="${ATDM_CONFIG_GNU_CXX_WARNINGS}"
  fi
  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}

elif [[ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0_OPENMPI-4.0.3" ]] ; then
  atdm_config_load_sparc_dev_module sparc-dev/gcc-7.2.0_openmpi-4.0.3
  unset OMP_NUM_THREADS  # SPARC module sets these and we must unset!
  unset OMP_PROC_BIND
  unset OMP_PLACES
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  if [[ "$ATDM_CONFIG_ENABLE_STRONG_WARNINGS" == "1" ]]; then
    export ATDM_CONFIG_CXX_FLAGS="${ATDM_CONFIG_GNU_CXX_WARNINGS}"
  fi
  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}
  export ATDM_CONFIG_MPI_EXEC=mpirun
  export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=-np
  export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

elif [ "$ATDM_CONFIG_COMPILER" == "INTEL-19.0.3_INTELMPI-2018.4" ]; then
  atdm_config_load_sparc_dev_module sparc-dev/intel-19.0.3_intelmpi-2018.4
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  if [[ "$ATDM_CONFIG_ENABLE_STRONG_WARNINGS" == "1" ]]; then
    export ATDM_CONFIG_CXX_FLAGS="${ATDM_CONFIG_INTEL_CXX_WARNINGS}"
  fi
  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}
  export ATDM_CONFIG_MPI_EXEC=mpirun
  export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=-np
  export ATDM_CONFIG_OPENMP_FORTRAN_FLAGS=-fopenmp
  export ATDM_CONFIG_OPENMP_FORTRAN_LIB_NAMES=gomp
  export ATDM_CONFIG_OPENMP_GOMP_LIBRARY=-lgomp
  # Point CMake 3.19 compiler checks to missing symbols
  export LDFLAGS="$LDFLAGS -lifcore"

elif [ "$ATDM_CONFIG_COMPILER" == "INTEL-19.0.3_MPICH2-3.2" ]; then
  atdm_config_load_sparc_dev_module sparc-dev/intel-19.0.3_mpich2-3.2
  export OMP_NUM_THREADS=3 # Because Si H. requested this
  export OMP_PROC_BIND=false
  unset OMP_PLACES
  export OMPI_CXX=`which icpc`
  export OMPI_CC=`which icc`
  export OMPI_FC=`which ifort`
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`
  if [[ "$ATDM_CONFIG_ENABLE_STRONG_WARNINGS" == "1" ]]; then
    export ATDM_CONFIG_CXX_FLAGS="${ATDM_CONFIG_INTEL_CXX_WARNINGS}"
  fi
  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}
  export ATDM_CONFIG_MPI_EXEC=mpirun
  export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=-np
  export ATDM_CONFIG_MPI_POST_FLAGS="-bind-to;none"
  export ATDM_CONFIG_OPENMP_FORTRAN_FLAGS=-fopenmp
  export ATDM_CONFIG_OPENMP_FORTRAN_LIB_NAMES=gomp
  export ATDM_CONFIG_OPENMP_GOMP_LIBRARY=-lgomp
  # Point CMake 3.19 compiler checks to missing symbols
  export LDFLAGS="$LDFLAGS -lifcore"

elif [ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1.243_GNU-7.2.0_OPENMPI-4.0.3" ]; then
  # ninja is running into issues with response files when building shared libraries with CUDA.
  # nvcc reports that no input files were given when generating shared libraries with nvcc.
  # Using the Unix Makefiles cmake generator works.
  export ATDM_CONFIG_USE_NINJA=OFF

  atdm_config_load_sparc_dev_module sparc-dev/cuda-10.1.243_gcc-7.2.0_openmpi-4.0.3

  # OpenMPI Settings
  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER}
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  # NOTE: The above export overrides the value set by the module load above
  export MPICC=`which mpicc`
  export MPICXX=`which mpicxx`
  export MPIF90=`which mpif90`

  # CUDA Settings
  if [[ ! -d /tmp/${USER} ]] ; then
    echo "Creating /tmp/${USER} for nvcc wrapper!"
    mkdir /tmp/${USER}
  fi

  export ATDM_CONFIG_MPI_EXEC=mpirun
  export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=-np
  export ATDM_CONFIG_MPI_POST_FLAGS="-bind-to;none"

  export ATDM_CONFIG_MKL_ROOT=${CBLAS_ROOT}
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return

fi

# ToDo: Update above to only load the compiler and MPI moudles and then
# directly set <TPL_NAME>_ROOT to point to the right TPLs.  This is needed to
# avoid having people depend on the SPARC modules.

# Set parallel test level based on OpenMP or not
if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
  #export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=$(($ATDM_CONFIG_MAX_NUM_CORES_TO_USE/2))
fi
# NOTE: Above, we use 1/2 as many parallel processes as cores on the machine
# to be safe.  Also, we need to set OMP_* env vars here because the SPARC
# modules change them!

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"
# NOTE: Above, we have to explicitly set the libs to use libbdf.so instead of
# libbdf.a because the former works and the latter does not and TriBITS is set
# up to only find static libs by default!

# BLAS and LAPACK

#export ATDM_CONFIG_BLAS_LIBS="-L${ATDM_CONFIG_MKL_ROOT}/mkl/lib/intel64;-L${ATDM_CONFIG_MKL_ROOT}/lib/intel64;-lmkl_intel_lp64;-lmkl_intel_thread;-lmkl_core;-liomp5"
#export ATDM_CONFIG_LAPACK_LIBS="-L${ATDM_CONFIG_MKL_ROOT}/mkl/lib/intel64"

# NOTE: The above does not work.  For some reason, the library 'iomp5' can't
# be found at runtime.  Instead, you have to explicitly list out the library
# files in order as shown below.  Very sad.

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${ATDM_CONFIG_MKL_ROOT}/mkl/lib/intel64 .so \
  mkl_intel_lp64 mkl_intel_thread mkl_core

atdm_config_add_libs_to_var ATDM_CONFIG_BLAS_LIBS ${ATDM_CONFIG_MKL_ROOT}/lib/intel64 .so \
  iomp5

export ATDM_CONFIG_LAPACK_LIBS=${ATDM_CONFIG_BLAS_LIBS}

# Boost

atdm_config_add_libs_to_var ATDM_CONFIG_BOOST_LIBS ${BOOST_ROOT}/lib .a \
  boost_program_options boost_system

# NOTE: Above, the SPARC-installed TPLs only have *.a files.  There are no
# *.so files.

# HDF5 and Netcdf

# Set the direct libs for HDF5 and NetCDF in case we use that option for
# building (see env var ATDM_CONFIG_USE_SPARC_TPL_FIND_SETTINGS).

if [[ "${ATDM_CONFIG_USE_MPI}" == "ON" ]] ; then
  USE_HDF5_ROOT="${HDF5_ROOT}"
  USE_NETCDF_ROOT="${NETCDF_ROOT}"
else
  USE_HDF5_ROOT="${SPARC_SERIAL_HDF5_ROOT}"
  USE_NETCDF_ROOT="${SPARC_SERIAL_NETCDF_ROOT}"
  unset HDF5_ROOT
  unset NETCDF_ROOT
  unset CGNS_ROOT
  # The SPARC TPL installer does not set the above ROOT vars for an non-MPI
  # (i.e. serial) build.  Therefore, we remove them here so that we can test
  # that the ATDMDevEnvSettings.cmake file can work correctly without these
  # set.
fi

export ATDM_CONFIG_HDF5_LIBS="-L${USE_HDF5_ROOT}/lib;${USE_HDF5_ROOT}/lib/libhdf5_hl.a;${USE_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

if [[ "${PNETCDF_ROOT}" == "" ]] ; then
  export PNETCDF_ROOT=${USE_NETCDF_ROOT}
fi
export ATDM_CONFIG_NETCDF_LIBS="-L${USE_NETCDF_ROOT}/lib;${USE_NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"

# SuperLUDist
if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
  # Set the default which is correct for all of the new TPL builds
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a
fi

# Finished!
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
