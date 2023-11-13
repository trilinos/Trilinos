################################################################################
#
# Set up env on weaver for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

# Handle compiler defaults

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA" ]] \
  || [[ "$ATDM_CONFIG_COMPILER" == "CUDA-11.2" ]] \
  || [[ "$ATDM_CONFIG_COMPILER" == "CUDA-11.2_GNU-8.3.0" ]] \
  || [[ "$ATDM_CONFIG_COMPILER" == "CUDA-11.2_GNU-8.3.0-OPENMPI-4.1.1" ]] \
   ; then
  export ATDM_CONFIG_COMPILER=CUDA-11.2-GNU-8.3.0-OPENMPI-4.1.1

else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  echo "*** Suppoted compilers include:"
  echo "***"
  echo "***   cuda-11.2-gnu-8.3.0  (cuda-11.2 default)"
  echo "***"
  return
fi

# Handle KOKKOS_ARCH

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] ; then
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

echo "Using weaver compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=OFF
export ATDM_CONFIG_USE_NINJA=OFF

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

source /etc/profile.d/modules.sh
# Temporary script to source while older modules are migrated to new system
source /projects/ppc64le-pwr9-rhel8/legacy-env.sh
module purge

module load git/2.10.1 python/3.7.3
# NOTE: Must load a git module since /usr/bin/git does not exist on the
# compute nodes.

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=20
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=40
fi

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]] ; then

  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-11.2-GNU-8.3.0-OPENMPI-4.1.1" ]] ; then
    module load cmake/3.24.2 cuda/11.2.2/gcc/8.3.1 openmpi/4.1.1/gcc/8.3.1/cuda/11.2.2 openblas/0.3.18/gcc/8.3.1 boost/1.70.0/gcc/8.3.1 metis/5.1.0/gcc/8.3.1 zlib/1.2.11/gcc/8.3.1 hdf5/1.10.7/gcc/8.3.1/openmpi/4.1.1 netcdf-c/4.8.1/gcc/8.3.1/openmpi/4.1.1 parallel-netcdf/1.12.2/gcc/8.3.1/openmpi/4.1.1 

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

  export ATDM_CONFIG_LAPACK_LIBS="-L${OPENBLAS_ROOT}/lib;-lopenblas"
  export ATDM_CONFIG_BLAS_LIBS="-L${OPENBLAS_ROOT}/lib;-lopenblas"
#  export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
#  export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"

  export ATDM_CONFIG_USE_CUDA=ON
#  export CUDA_LAUNCH_BLOCKING=1
#  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
#  export KOKKOS_NUM_DEVICES=2
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=2
  # Avoids timeouts due to not running on separate GPUs (e.g. see #2446)

fi

# CMake
module load cmake/3.24.2

# HWLOC

export ATDM_CONFIG_USE_HWLOC=OFF

# Let's see if the TPLs loaded by devpack/20180517/openmpi/2.1.2/gcc/7.2.0/cuda/9.2.88 work for SPARC?
#export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.a;${BINUTILS_ROOT}/lib/libiberty.a"

# HDF5 and Netcdf

# NOTE: HDF5_ROOT and NETCDF_C_ROOT should already be set in env from above
# module loads!

# However, set the direct libs for HDF5 and NetCDF in case we use that option
# for building (see env var ATDM_CONFIG_USE_SPARC_TPL_FIND_SETTINGS).

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

if [[ "${PARALLEL_NETCDF_ROOT}" == "" ]] ; then
  export NETCDF_ROOT=${NETCDF_C_ROOT}
  export PNETCDF_ROOT=${NETCDF_C_ROOT}
  export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_C_ROOT}/lib;${NETCDF_C_ROOT}/lib/libnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"
else
  export NETCDF_ROOT=${NETCDF_C_ROOT}
  export PNETCDF_ROOT=${PARALLEL_NETCDF_ROOT}
  export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_C_ROOT}/lib;${NETCDF_C_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"
fi


# SuperLUDist

#if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
#  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
#  export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib/libsuperlu_dist.a;${METIS_ROOT}/lib/libmetis.a"
#fi

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
