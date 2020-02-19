################################################################################
#
# Set up env on ats2 for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

# ats2 jobs all use the same environmnet changes to the
# sourced script below will impact jobs on both of those
# machines. please be mindful of this when making changes

#
# Handle KOKKOS_ARCH
#

if [[ "$ATDM_CONFIG_COMPILER" == "GNU"* || \
  "$ATDM_CONFIG_COMPILER" == "XL"* ]]; then
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9
    sparc_tpl_arch=pwr9
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
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Volta70" ]] ; then
    export ATDM_CONFIG_KOKKOS_ARCH=Power9,Volta70
    sparc_tpl_arch=v100
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

echo "Using $ATDM_CONFIG_SYSTEM_NAME compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

# Some basic settings
export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=8
# ToDo: Experiment with increasing the build parallelism for different
# compilers, different build configurations, on the login node vs. the compute
# node, etc.

# Set ctest -j parallel level for non-CUDA builds
if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

# NOTE: We do *NOT* purge the modules first like on other systems because that
# messes up the ATS-2 env.  Therefore, it is recommended that the user load a
# new shell and then load one of these envs.

# Load common modules for all builds
module load git/2.20.0
module load cmake/3.14.5

#
# Load compiler modules, TPL modules, and point to SPARC TPL install base dirs
#

if [[ "$ATDM_CONFIG_COMPILER" == *"GNU-7.3.1_SPMPI-2019.06.24" ]]; then
  module load gcc/7.3.1
  module load lapack/3.8.0-gcc-4.9.3

  export CBLAS_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-gcc-4.9.3
  export LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-gcc-4.9.3
  export COMPILER_ROOT=/usr/tce/packages/gcc/gcc-7.3.1
  export SPARC_HDF5=hdf5-1.10.5

  # eharvey: TODO: remove COMPILER_ROOT and other unused exports below.
  export PATH=${COMPILER_ROOT}/bin:${PATH}
  export LD_LIBRARY_PATH=${COMPILER_ROOT}/lib:${LD_LIBRARY_PATH}
  export BINUTILS_ROOT=${COMPILER_ROOT}
  export LIBRARY_PATH=${BINUTILS_ROOT}/lib
  export LIBRARY_PATH=${CBLAS_ROOT}/lib:${LIBRARY_PATH}
  export INCLUDE=${BINUTILS_ROOT}/include:${INCLUDE}
  export CPATH=${BINUTILS_ROOT}/include:${CPATH}

  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1.243_"* ]]; then
    sparc_tpl_ext=ats2-${sparc_tpl_arch}_cuda-10.1.243_gcc-7.3.1
    sparc_tpl_mpi_ext=ats2-${sparc_tpl_arch}_cuda-10.1.243_gcc-7.3.1_spmpi-2019.06.24
  else
    sparc_tpl_ext=ats2-${sparc_tpl_arch}_gcc-7.3.1
    sparc_tpl_mpi_ext=ats2-${sparc_tpl_arch}_gcc-7.3.1_spmpi-2019.06.24
  fi

elif [[ "$ATDM_CONFIG_COMPILER" == *"XL-2019.08.20_SPMPI-2019.06.24_DISABLED" ]]; then
  module load xl/2019.08.20
  module load lapack/3.8.0-xl-2019.08.20
  module load gmake/4.2.1

  # Ninja not available for XL until cmake 3.17.0
  export ATDM_CONFIG_USE_NINJA=OFF

  export CBLAS_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-P9-xl-2019.08.20
  export LAPACK_ROOT=/usr/tcetmp/packages/lapack/lapack-3.8.0-P9-xl-2019.08.20
  export COMPILER_ROOT=/usr/tce/packages/xl/xl-2019.08.20
  export SPARC_HDF5=hdf5-1.8.20

  # eharvey: TODO: remove COMPILER_ROOT and other exports below.
  export PATH=${COMPILER_ROOT}/bin:${PATH}
  export LD_LIBRARY_PATH=${COMPILER_ROOT}/lib:${LD_LIBRARY_PATH}
  export BINUTILS_ROOT=/usr/tce/packages/gcc/gcc-7.3.1
  export LIBRARY_PATH=${BINUTILS_ROOT}/lib
  export LIBRARY_PATH=${CBLAS_ROOT}/lib:${LIBRARY_PATH}
  export INCLUDE=${BINUTILS_ROOT}/include:${INCLUDE}
  export CPATH=${BINUTILS_ROOT}/include:${CPATH}

  if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1.243_"* ]]; then
    export LD_LIBRARY_PATH=${BINUTILS_ROOT}/rh/lib/gcc/ppc64le-redhat-linux/7:${LD_LIBRARY_PATH}
    sparc_tpl_ext=ats2-${sparc_tpl_arch}_cuda-10.1.243_xl-2019.08.20
    sparc_tpl_mpi_ext=ats2-${sparc_tpl_arch}_cuda-10.1.243_xl-2019.08.20_spmpi-2019.06.24
  else
    sparc_tpl_ext=ats2-${sparc_tpl_arch}_xl-2019.08.20
    sparc_tpl_mpi_ext=ats2-${sparc_tpl_arch}_xl-2019.08.20_spmpi-2019.06.24
  fi

else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
  echo "***"
  return

fi

#
# Load module and do other setup for CUDA bulids
#

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA-10.1.243_"* ]]; then

  module load cuda/10.1.243

  # OpenMPI Settings
  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER}
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  export OMPI_CC=`which gcc`
  export OMPI_FC=`which gfortran`
  export LLNL_USE_OMPI_VARS="y"

  # CUDA Settings
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
  export TMPDIR=/tmp/$(whoami)
  if [[ ! -d ${TMPDIR} ]] ; then
    echo "Creating ${TMPDIR} for nvcc wrapper!"
    mkdir ${TMPDIR}
  fi
  # ATDM Settings
  export ATDM_CONFIG_USE_CUDA=ON
  export ATDM_CONFIG_USE_OPENMP=OFF
  export ATDM_CONFIG_USE_PTHREADS=OFF
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=4

  # Kokkos Settings
  export ATDM_CONFIG_Kokkos_ENABLE_SERIAL=ON
  export KOKKOS_NUM_DEVICES=4

  # CTEST Settings
  # TPETRA_ASSUME_CUDA_AWARE_MPI is used by cmake/std/atdm/ats2/trilinos_jsrun
  export TPETRA_ASSUME_CUDA_AWARE_MPI=0
  # Trilinos_CTEST_RUN_CUDA_AWARE_MPI is used by cmake/ctest/driver/atdm/ats2/local-driver.sh
  export Trilinos_CTEST_RUN_CUDA_AWARE_MPI=1

elif [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then

  echo
  echo "***"
  echo "*** ERROR: CUDA version in COMPILER=$ATDM_CONFIG_COMPILER"
  echo "*** is not supported on this system!  Only CUDA-10.1.243"
  echo "*** is currently supported!"
  echo "***"
  return

fi

#
# Final setup for all build configurations
#

# Common module - requires compiler to be loaded first
module load spectrum-mpi/2019.06.24

# Prepend path to ninja after all of the modules are loaded
export PATH=/projects/atdm_devops/vortex/ninja-fortran-1.8.2:$PATH

# ATDM specific config variables
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"

# NOTE: Invalid libbfd.so requires below for Trilinos to compile
export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.a;-lz;${BINUTILS_ROOT}/lib/libiberty.a"

export ATDM_CONFIG_SPARC_TPL_BASE=/projects/sparc/tpls/ats2-${sparc_tpl_arch}

sparc_tpl_base=${ATDM_CONFIG_SPARC_TPL_BASE}

# Commont ROOT config variables
export BOOST_ROOT=${sparc_tpl_base}/boost-1.65.1/00000000/${sparc_tpl_ext}
export HDF5_ROOT=${sparc_tpl_base}/hdf5-1.10.5/00000000/${sparc_tpl_mpi_ext}
export CGNS_ROOT=${sparc_tpl_base}/cgns-c09a5cd/27e5681f1b74c679b5dcb337ac71036d16c47977/${sparc_tpl_mpi_ext}
export PNETCDF_ROOT=${sparc_tpl_base}/pnetcdf-1.10.0/6144dc67b2041e4093063a04e89fc1e33398bd09/${sparc_tpl_mpi_ext}
export NETCDF_ROOT=${sparc_tpl_base}/netcdf-4.7.0/58bc48d95be2cc9272a18488fea52e1be1f0b42a/${sparc_tpl_mpi_ext}
export PARMETIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export METIS_ROOT=${sparc_tpl_base}/parmetis-4.0.3/00000000/${sparc_tpl_mpi_ext}
export SUPERLUDIST_ROOT=${sparc_tpl_base}/superlu_dist-5.4.0/a3121eaff44f7bf7d44e625c3b3d2a9911e58876/${sparc_tpl_mpi_ext}

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${SEMS_PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

if [[ "${ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS}" == "" ]] ; then
  export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
  export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a"
fi

# Set common MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_EXEC=${ATDM_SCRIPT_DIR}/ats2/trilinos_jsrun

export ATDM_CONFIG_MPI_POST_FLAGS="--rs_per_socket;4"
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="-p"

# Set common default compilers
export CC=mpicc
export CXX=mpicxx
export F77=mpifort
export FC=mpifort
export F90=mpifort

# Define function atdm_run_script_on_compute_node
unset atdm_run_script_on_compute_node

source $ATDM_SCRIPT_DIR/common/define_run_on_lsf_compute_node_func.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE

export ATDM_CONFIG_TRIL_CMAKE_INSTALL_PREFIX_DATE_BASE_DEFAULT=/projects/atdm_devops/trilinos_installs/
