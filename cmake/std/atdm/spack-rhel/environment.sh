################################################################################
#
# Set up env using spack build tpls on Linux RHEL platforms for ATMD builds of
# Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

#
# Functions
#

function load_spack_tpl_modules() {
  compiler_dash_version=$1
  mpi_dash_version=$2

  module load spack-netlib-lapack/3.8.0-${compiler_dash_version}
  export BLAS_ROOT=$NETLIB_LAPACK_ROOT
  export LAPACK_ROOT=$NETLIB_LAPACK_ROOT

  module load spack-binutils/2.31.1-${compiler_dash_version}
  module load spack-boost/1.59.0-${compiler_dash_version}
  module load spack-superlu/4.3-${compiler_dash_version}
  module load spack-zlib/1.2.11-${compiler_dash_version}
  module load spack-metis/5.1.0-${compiler_dash_version}
  module load spack-hdf5/1.8.21-${compiler_dash_version}-${mpi_dash_version}
  module load spack-netcdf/4.4.1-${compiler_dash_version}-${mpi_dash_version}
  module load spack-parmetis/4.0.3-${compiler_dash_version}-${mpi_dash_version}
  module load spack-cgns/snl-atdm-${compiler_dash_version}-${mpi_dash_version}
  module load spack-superlu-dist/5.4.0-${compiler_dash_version}-${mpi_dash_version}

  # NOTE: Above, we don't need to explicilty load the module for
  # sparc-parallel-netcdf because it is implicitly enabled by loading the
  # module for spack-netcdf.  And we not want to explicitly load any modules
  # for packages that we have not set the explicit versions on as that breaks
  # when we upgarde spack (and spack uses a new version).

}


#
# Allow KOKKOS_ARCH to be set but unset it if DEFAULT
#

if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  unset ATDM_CONFIG_KOKKOS_ARCH
else
  echo
  echo "***"
  echo "*** ERROR: Specifying KOKKOS_ARCH is not supported on RHEL ATDM builds"
  echo "*** remove '$ATDM_CONFIG_KOKKOS_ARCH' from ATDM_CONFIG_BUILD_NAME=$ATDM_CONFIG_BUILD_NAME"
  echo "***"
  return
fi

echo "Using SPACK RHEL compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
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
module load spack-git/2.20.1
module load spack-cmake/3.14.5
module load spack-ninja-fortran/1.9.0.2.g99df1
module load spack-python/2.7.15  # For EMPIRE tests

# Assume binutils lib dir is 'lib64' by default
BINUTILS_LIBIBERTY_LIB_DIR_NAME=lib64

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

if [ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0_OPENMPI-1.10.1" ]; then

  module load spack-gcc/7.2.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`

  module load spack-openmpi/1.10.1-gcc-7.2.0

  load_spack_tpl_modules gcc-7.2.0 openmpi-1.10.1

elif [ "$ATDM_CONFIG_COMPILER" == "GNU-7.2.0_OPENMPI-2.1.2" ]; then

  module load spack-gcc/7.2.0
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`

  module load spack-openmpi/2.1.2-gcc-7.2.0

  load_spack_tpl_modules gcc-7.2.0 openmpi-2.1.2

elif [ "$ATDM_CONFIG_COMPILER" == "CLANG-5.0.1_OPENMPI-1.10.2" ]; then

  module load spack-gcc/4.9.3  # Need to find gfortran and for llvm to find the right gcc!
  module load spack-llvm/5.0.1
  export OMPI_CXX=`which clang++`
  export OMPI_CC=`which clang`
  export OMPI_FC=`which gfortran`

  module load spack-openmpi/1.10.2-gcc-4.9.3

  load_spack_tpl_modules clang-5.0.1 openmpi-1.10.2

  BINUTILS_LIBIBERTY_LIB_DIR_NAME=lib

  export ATDM_CONFIG_EXTRA_LINK_FLAGS="${GCC_ROOT}/lib64/libgfortran.so"
  # NOTE: Above, for some reason, need to provide full path to libgfortran as
  # -L and -l don't work.

  # This does NOT work so we have to use above instead!
  #export ATDM_CONFIG_EXTRA_LINK_FLAGS="-L${GCC_ROOT}/lib64;-lgfortran"

else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return
fi

if [[ "${ATDM_CONFIG_SHARED_LIBS}" == "ON" ]] ; then
  ATDM_CONFIG_TPL_LIB_EXT=so
else
  ATDM_CONFIG_TPL_LIB_EXT=a
fi

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=-lhwloc

export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib64;-llapack"
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib64;-lblas"

export ATDM_CONFIG_HDF5_LIBS="${HDF5_ROOT}/lib/libhdf5_hl.so;${HDF5_ROOT}/lib/libhdf5.so;${ZLIB_ROOT}/lib/libz.so;-ldl"
#echo ATDM_CONFIG_HDF5_LIBS=$ATDM_CONFIG_HDF5_LIBS

export PNETCDF_ROOT=${SEMS_NETCDF_ROOT}

export ATDM_CONFIG_METIS_LIBS="${METIS_ROOT}/lib/libmetis.so"
export ATDM_CONFIG_PARMETIS_LIBS="${METIS_ROOT}/lib/libmetis.so;${PARMETIS_ROOT}/lib/libparmetis.so"
#export ATDM_CONFIG_PNETCDF_LIBS=
export ATDM_CONFIG_CGNS_LIBS="${CGNS_ROOT}/lib/libcgns.so;${ATDM_CONFIG_HDF5_LIBS}"

#export METIS_LIBRARY_DIRS=${METIS_ROOT}/lib
#export METIS_INCLUDE_DIRS=${METIS_ROOT}/include

export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"
# NOTE: SEMS does not provide the correct *.so files for NetCDF so we can't
# use them in a shared lib build :-(

export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.so;${BINUTILS_ROOT}/${BINUTILS_LIBIBERTY_LIB_DIR_NAME}/libiberty.a;${GETTEXT_ROOT}/lib/libintl.so;${LIBICONV_ROOT}/lib/libiconv.so"
# NOTE: Above, we have to explicitly set the libs to use libbdf.so instead of
# libbdf.a because the former works and the latter does not and TriBITS is set
# up to only find static libs by default!
unset BINUTILS_ROOT
# NOTE: Above, you have to unset BINUTILS_ROOT or the CMake configure of
# Trilinos dies in the compiler check.  (We should investigate this more at
# some point.)

# SuperLUDist
if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
  # Set the default which is correct for all of the new TPL builds
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLU_DIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLU_DIST_ROOT}/lib/libsuperlu_dist.so
fi


# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
