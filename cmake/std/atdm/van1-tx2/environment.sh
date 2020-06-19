################################################################################
#
# Set up env on a ARM ATSE builds of Trilinos
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

echo "Using ARM ATSE compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

#
# Assert and set KOKKOS_ARCH
#

if   [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ]] \
  || [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "TX2" ]] \
  || [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "" ]] \
  ; then
  export ATDM_CONFIG_KOKKOS_ARCH=ARMv8-TX2
else
  echo
  echo "***"
  echo "*** ERROR: KOKKOS_ARCH='${ATDM_CONFIG_KOKKOS_ARCH}' was parsed from the the buildname '${ATDM_CONFIG_BUILD_NAME}'.  Only one KOKKOS_ARCH is supported for this system.  Please remove that KOKKOS_ARCH keyword from the buildname!"
  echo "***"
  return
fi

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

export ATDM_CONFIG_BUILD_COUNT=40  # Assume building on a compute node!
if [[ "${ATDM_CONFIG_BUILD_TYPE}" == "DEBUG" ]] ; then
  export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=20
  # Above: The 'dbg' build on 'stria' is randomly failing the link of some ROL
  # execuables due to running out of memory when using 40 parallel link jobs.
  # Reducing this is to avoid that.  See CDOFA-117.
fi

#
# Load the modules
#

module purge

if [[ "$ATDM_CONFIG_COMPILER" == "ARM-20.0_OPENMPI-4.0.2" ]]; then
  module load devpack-arm
  module unload yaml-cpp
  # provides numpy module for empire
  module load python/3.6.8
  module load arm/20.0
  # Check if openmpi is already loaded. If it is, swap it. Otherwise, just load ompi4.
  if [[ ! -z $(module list openmpi | grep '1)' | awk -F ' ' '{print $2}') ]]; then
    module swap $(module list openmpi | grep '1)' | awk -F ' ' '{print $2}') openmpi4/4.0.2
  else
    module load openmpi4/4.0.2
  fi
  module load armpl/20.0.0

  export LAPACK_ROOT="$ARMPL_LIB"
  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-larmpl_ilp64_mp"
  export ATDM_CONFIG_BLAS_LIBS="-L${LAPACK_ROOT};-larmpl_ilp64_mp"

  # We'll use TPL_ROOT for consistency across ATDM environments
  export MPI_ROOT=${MPI_DIR}
  export BLAS_ROOT=${ARMPL_DIR}
  export LAPACK_ROOT=${ARMPL_DIR}
  export HDF5_ROOT=${HDF5_DIR}
  export NETCDF_ROOT=${NETCDF_DIR}
  export PNETCDF_ROOT=${PNETCDF_DIR}
  export ZLIB_ROOT=${ZLIB_DIR}
  export CGNS_ROOT=${CGNS_DIR}
  export BOOST_ROOT=${BOOST_DIR}
  export METIS_ROOT=${METIS_DIR}
  export PARMETIS_ROOT=${PARMETIS_DIR}
  export SUPERLUDIST_ROOT=${SUPERLU_DIST_DIR}
  export BINUTILS_ROOT=${BINUTILS_DIR}

  module load git/2.19.2
elif [[ "$ATDM_CONFIG_COMPILER" == "ARM-20.1_OPENMPI-4.0.3" ]]; then
  module load sparc-dev/arm-20.1_openmpi-4.0.3

  if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
    unset OMP_PLACES
    unset OMP_PROC_BIND
  fi

  # We'll use TPL_ROOT for consistency across ATDM environments
  export MPI_ROOT=${MPI_DIR}
  export BLAS_ROOT=${ARMPL_DIR}
  export HDF5_ROOT=${HDF5_DIR}
  export NETCDF_ROOT=${NETCDF_DIR}
  export PNETCDF_ROOT=${PNETCDF_DIR}
  export ZLIB_ROOT=${ZLIB_DIR}
  export CGNS_ROOT=${CGNS_DIR}
  export METIS_ROOT=${METIS_DIR}
  export PARMETIS_ROOT=${PARMETIS_DIR}
  export SUPERLUDIST_ROOT=${SUPERLU_DIST_DIR}
  export BINUTILS_ROOT=${BINUTILS_DIR}
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

if [[ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ]] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=32
  export OMP_PROC_BIND=FALSE
  export OMP_NUM_THREADS=1
fi

# Common modules for all builds
module load ninja
module load cmake/3.12.2

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=

# BINUTILS settings
export ATDM_CONFIG_BINUTILS_LIBS="-L${BINUTILS_ROOT}/lib;-lbfd"

# CGNS settings
export ATDM_CONFIG_CGNS_LIBRARY_NAMES="cgns"

# HDF5 settings
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

# NETCDF settings
export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"

# BLAS settings
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-larmpl_lp64_mp;-larmflang;-lomp"

# LAPACK settings
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-larmpl_lp64_mp;-larmflang;-lomp"

# SuperLUDist settings
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS="${SUPERLUDIST_ROOT}/include"
export ATDM_CONFIG_SUPERLUDIST_LIBS="-L${SUPERLUDIST_ROOT}/lib;-lsuperlu_dist"

# METIS settings - force linker order!
export ATDM_CONFIG_PARMETIS_LIBS="${PARMETIS_ROOT}/lib/libparmetis.a;${METIS_ROOT}/lib/libmetis.a"

#
# MPI Settings
#

# Set common MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

# Set common default compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F90=mpif90

# MPI_EXEC settings
export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="-np"
export ATDM_CONFIG_MPI_EXEC="mpirun"
export ATMD_CONFIG_MPI_USE_COMPILER_WRAPPERS=ON

export ATDM_CONFIG_WCID_ACCOUNT_DEFAULT=fy150090

#
# Done
#
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
