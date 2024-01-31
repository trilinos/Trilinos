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

if [[ "${ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE}" != "1" ]] ; then
  module purge
else
  echo "NOTE: ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE=1 is set so using pre-loaded sparc-dev module!"
fi

if [[ "$ATDM_CONFIG_COMPILER" == "ARM-20.1_OPENMPI-4.0.5" ]]; then
  module load sparc-dev/arm-20.1_openmpi-4.0.5
  # devpack includes cmake/3.17.1 which is too old; swap it for this version
  module load cmake/3.27.4
  module unload yaml-cpp
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
  # SUPERLU_DIST and NETCDF libraries live in $ROOT/lib in this version
  export ATDM_CONFIG_SUPERLUDIST_LIBS="-L${SUPERLUDIST_ROOT}/lib;-lsuperlu_dist"
  export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"
elif [[ "$ATDM_CONFIG_COMPILER" == "ARM-22.1_OPENMPI-4.0.5" ]]; then
  export BLAS_ROOT=${ARMPL_DIR}
  module load sparc-dev/arm-22.1_openmpi-4.0.5
  module unload sparc-cmake/3.23.2
  module load cmake/3.27.4
  module unload yaml-cpp
  # These TPLs are not included in the devpack:
  module load superlu/5.2.1
  module load superlu_dist/5.4.0
  module load binutils/2.33.1
  module load boost/1.72.0
  module load phdf5/1.10.5
  module load netcdf/4.6.3
  module load pnetcdf/1.11.1
  module load metis/5.1.0
  module load parmetis/4.0.3
  module load cgns/3.4.0
  # Most <TPL>_ROOT paths are not set by TPLs, instead we have
  # <TPL>_DIR. Trilinos looks for the _ROOT form.
  export MPI_ROOT=${MPI_DIR}
  export HDF5_ROOT=${HDF5_DIR}
  export NETCDF_ROOT=${NETCDF_DIR}
  export PNetCDF_ROOT=${PNETCDF_DIR}
  export PNETCDF_ROOT=${PNETCDF_DIR}
  export SUPERLU_ROOT=${SUPERLU_DIR}
  export SUPERLUDist_ROOT=${SUPERLU_DIST_DIR}
  export SuperLUDist_ROOT=${SUPERLU_DIST_DIR}
  export METIS_ROOT=${METIS_DIR}
  export PARMETIS_ROOT=${PARMETIS_DIR}
  export CGNS_ROOT=${CGNS_DIR}
  export BINUTILS_ROOT=${BINUTILS_DIR}
  export ATDM_CONFIG_SUPERLUDIST_LIBS="-L${SUPERLUDIST_ROOT}/lib;-lsuperlu_dist"
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported!"
  echo "***"
  return
fi

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  unset OMP_PLACES
  unset OMP_PROC_BIND
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

export ATDM_CONFIG_USE_HWLOC=OFF
export HWLOC_LIBS=

# BINUTILS settings
export ATDM_CONFIG_BINUTILS_LIBS="-L${BINUTILS_ROOT}/lib;-lbfd"

# CGNS settings
export ATDM_CONFIG_CGNS_LIBRARY_NAMES="cgns"

# HDF5 settings
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

# BLAS settings
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-larmpl_lp64_mp;-larmflang;-lomp"

# LAPACK settings
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-larmpl_lp64_mp;-larmflang;-lomp"

# Common SuperLUDist settings
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS="${SUPERLUDIST_ROOT}/include"

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
export ATDM_CONFIG_MPI_EXEC="$MPI_ROOT/bin/mpirun"
export ATMD_CONFIG_MPI_USE_COMPILER_WRAPPERS=ON

export ATDM_CONFIG_WCID_ACCOUNT_DEFAULT=fy150090

#
# Done
#
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
