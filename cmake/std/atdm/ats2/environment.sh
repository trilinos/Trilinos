################################################################################
#
# Set up env on ATS2 systems for supported builds
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
  module load spack-gettext/0.19.8.1-${compiler_dash_version} # for binutils
  module load spack-libiconv/1.15-${compiler_dash_version} # for gettext 
  module load spack-boost/1.59.0-${compiler_dash_version}
  module load spack-superlu/4.3-${compiler_dash_version}
  module load spack-zlib/1.2.11-${compiler_dash_version}
  module load spack-metis/5.1.0-${compiler_dash_version}
  module load spack-hdf5/1.8.21-${compiler_dash_version}-${mpi_dash_version}
  module load spack-netcdf/4.4.1-${compiler_dash_version}-${mpi_dash_version}
  module load spack-parallel-netcdf/1.11.0-${compiler_dash_version}-${mpi_dash_version}
  module load spack-parmetis/4.0.3-${compiler_dash_version}-${mpi_dash_version}
  module load spack-cgns/snl-atdm-${compiler_dash_version}-${mpi_dash_version}
  module load spack-superlu-dist/6.1.0-${compiler_dash_version}-${mpi_dash_version}

}


# Don't allow KOKKOS_ARCH to be set
if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  unset ATDM_CONFIG_KOKKOS_ARCH
else
  echo
  echo "***"
  echo "*** ERROR: Specifying KOKKOS_ARCH is not supported on the 'ats2' env!"
  echo "*** Remove '$ATDM_CONFIG_KOKKOS_ARCH' from ATDM_CONFIG_BUILD_NAME=$ATDM_CONFIG_BUILD_NAME"
  echo "***"
  return
fi

echo "Using ATS-2 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

# Set the build and link paralel job levels
if   [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] \
  && [[ "${ATDM_CONFIG_CUDA_RDC}" == "ON" ]] \
  ; then
  if [[ "${ATDM_CONFIG_BUILD_TYPE}" == *"DEBUG" ]] ;then
    export ATDM_CONFIG_BUILD_COUNT=20
    export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=10
  else
    export ATDM_CONFIG_BUILD_COUNT=32
    export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=16
  fi
  # When CUDA+RDC is enabled, using all 64 cores to build and link results in
  # build errors as described in #4502.
else
  export ATDM_CONFIG_BUILD_COUNT=64
fi

# Set the parallel testing level
if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=20
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=40
fi

# Load base modules modules 
module purge
module load spack-cmake/3.13.4-gcc-7.3.1  # ToDo: Remove compiler from name of these!
module load spack-git/2.20.1-gcc-7.3.1
module load spack-ninja-fortran/1.7.2.gaad58-gcc-7.3.1

if [ "$ATDM_CONFIG_COMPILER" == "GNU-7.3.1_OPENMPI-2.1.2" ]; then

  module load gcc/7.3.1
  export OMPI_CXX=`which g++`
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`

  module load spack-openmpi/2.1.2-gcc-7.3.1

  load_spack_tpl_modules gcc-7.3.1 openmpi-2.1.2

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

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;${ZLIB_ROOT}/lib/libz.a;-ldl"
#echo ATDM_CONFIG_HDF5_LIBS=$ATDM_CONFIG_HDF5_LIBS

export PNETCDF_ROOT=${PARALLEL_NETCDF_ROOT}

export ATDM_CONFIG_METIS_LIBS="${METIS_ROOT}/lib/libmetis.so"
export ATDM_CONFIG_PARMETIS_LIBS="${METIS_ROOT}/lib/libmetis.so;${PARMETIS_ROOT}/lib/libparmetis.so"
#export ATDM_CONFIG_PNETCDF_LIBS=
export ATDM_CONFIG_CGNS_LIBS="${CGNS_ROOT}/lib/libcgns.a;-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

#export METIS_LIBRARY_DIRS=${METIS_ROOT}/lib
#export METIS_INCLUDE_DIRS=${METIS_ROOT}/include

export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.${ATDM_CONFIG_TPL_LIB_EXT};${BOOST_ROOT}/lib/libboost_system.${ATDM_CONFIG_TPL_LIB_EXT};${NETCDF_ROOT}/lib/libnetcdf.a;${PARALLEL_NETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"
# NOTE: SEMS does not provide the correct *.so files for NetCDF so we can't
# use them in a shared lib build :-(

export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.so;${BINUTILS_ROOT}/lib64/libiberty.a;${GETTEXT_ROOT}/lib/libintl.a;${LIBICONV_ROOT}/lib/libiconv.so"
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

export ATDM_CONFIG_MPI_POST_FLAGS="-map-by;socket:PE=4"

#
# Done!
#

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
