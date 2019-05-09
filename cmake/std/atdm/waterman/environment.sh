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

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] && \
  [[ "${ATDM_CONFIG_CUDA_RDC}" == "ON" ]] ; then
  export ATDM_CONFIG_BUILD_COUNT=32
  export ATDM_CONFIG_PARALLEL_LINK_JOBS_LIMIT=16
  # When CUDA+RDC is enabled, using all 64 cores to build and link results in
  # build errors as described in #4502.
else
  export ATDM_CONFIG_BUILD_COUNT=64
fi

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

# CMake and ninja
module swap cmake/3.6.2 cmake/3.12.3
module load ninja/1.7.2

# HWLOC

export ATDM_CONFIG_USE_HWLOC=OFF

# Let's see if the TPLs loaded by devpack/20180517/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 work for SPARC?

export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.a;${BINUTILS_ROOT}/lib/libiberty.a"

#CGNS_ROOT=/home/projects/sparc/tpls/waterman/cgns-develop/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#HDF5_ROOT=/home/projects/sparc/tpls/waterman/hdf5-1.8.20/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#METIS_ROOT=/home/projects/sparc/tpls/waterman/parmetis-4.0.3/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#NETCDF_ROOT=/home/projects/sparc/tpls/waterman/netcdf-4.6.1/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#PARMETIS_ROOT=/home/projects/sparc/tpls/waterman/parmetis-4.0.3/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#PNETCDF_ROOT=/home/projects/sparc/tpls/waterman/pnetcdf-1.10.0/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#SGM_ROOT=/home/projects/sparc/tpls/waterman/sgm-develop/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2
#SUPERLUDIST_ROOT=/home/projects/sparc/tpls/waterman/superlu_dist-4.2/waterman-gpu_gcc-7.2.0_cuda-9.2.88_openmpi-2.1.2

# HDF5 and Netcdf

# NOTE: HDF5_ROOT and NETCDF_ROOT should already be set in env from above
# module loads!

# However, set the direct libs for HDF5 and NetCDF in case we use that option
# for building (see env var ATDM_CONFIG_USE_SPARC_TPL_FIND_SETTINGS).

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${SEMS_PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

#export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
#export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

# SuperLUDist

if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib/libsuperlu_dist.a;${METIS_ROOT}/lib/libmetis.a"
fi

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAGS="-map-by;socket:PE=4"

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
