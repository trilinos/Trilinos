################################################################################
#
# Set up env on ats2 for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

# ats2 jobs all use the same environment changes to the
# sourced script below will impact jobs on both of those
# machines. please be mindful of this when making changes

#
# Handle KOKKOS_ARCH
#

if [[ "$ATDM_CONFIG_COMPILER" == "GNU"* ]]; then
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" ]] ; then
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
  if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Power9" || \
    "$ATDM_CONFIG_KOKKOS_ARCH" == "Volta70" ]] ; then
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

echo "Using $ATDM_CONFIG_SYSTEM_NAME compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

# Some basic settings
export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_BUILD_COUNT=8 # Assume building on the shared login node!

# Set ctest -j parallel level for non-CUDA builds
if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

# Purge then load StdEnv to get back to a fresh env in case previous other
# modules were loaded.
if [[ "${ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE}" != "1" ]] ; then
  module purge --silent
  module load StdEnv
else
  echo "NOTE: ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE=1 is set so using pre-loaded sparc-dev module!"
fi

# Load the sparc-dev/xxx module
sparc_module_name=$(get_sparc_dev_module_name "$ATDM_CONFIG_COMPILER")
atdm_config_load_sparc_dev_module ${sparc_module_name}

# Set up stuff related the the host compiler

if [[ "$ATDM_CONFIG_COMPILER" == *"GNU"* ]]; then

  export COMPILER_ROOT=/usr/tce/packages/gcc/gcc-7.3.1
  export LD_LIBRARY_PATH=${COMPILER_ROOT}/lib:${LD_LIBRARY_PATH}
  export BINUTILS_ROOT=${COMPILER_ROOT}
  export LIBRARY_PATH=${BINUTILS_ROOT}/lib
  export LIBRARY_PATH=${CBLAS_ROOT}/lib:${LIBRARY_PATH}
  export INCLUDE=${BINUTILS_ROOT}/include:${INCLUDE}
  export CPATH=${BINUTILS_ROOT}/include:${CPATH}

  export ATDM_CONFIG_USE_NINJA=ON

fi

# Set up stuff related to CUDA

if [[ "$ATDM_CONFIG_COMPILER" == "CUDA"* ]]; then

  export CUDA_BIN_PATH=$CUDA_HOME

  # OpenMPI Settings
  export OMPI_CXX=${ATDM_CONFIG_NVCC_WRAPPER}
  if [ ! -x "$OMPI_CXX" ]; then
      echo "No nvcc_wrapper found"
      return
  fi
  # NOTE: The above export overrides the value set by the module load above

  # CUDA Settings
  if [[ ! -d /tmp/${USER} ]] ; then
    echo "Creating /tmp/${USER} for nvcc wrapper!"
    mkdir /tmp/${USER}
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
  # Trilinos_CTEST_RUN_GPU_AWARE_MPI is used by ats2/local-driver.sh
  export Trilinos_CTEST_RUN_GPU_AWARE_MPI=1

fi

#
# Final setup for all build configurations
#

# Prepend path to ninja after all of the modules are loaded
export PATH=/projects/atdm_devops/vortex/ninja-fortran-1.8.2:$PATH

# Set a standard git so everyone has the same git
module load git/2.20.0

# ATDM specific config variables
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
export ATDM_CONFIG_BLAS_LIBS="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"

# NOTE: Invalid libbfd.so requires below for Trilinos to compile
export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib/libbfd.a;-lz;${BINUTILS_ROOT}/lib/libiberty.a"

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib64;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"

export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a"

# Set common MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_EXEC=${ATDM_SCRIPT_DIR}/ats2/trilinos_jsrun

export ATDM_CONFIG_MPI_POST_FLAGS="--rs_per_socket;4"
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="-p"

if [[ "${ATDM_CONFIG_COMPLEX}" == "ON" ]] ; then
  export ATDM_CONFIG_MPI_PRE_FLAGS="-M;-mca coll ^ibm"
  # NOTE: We have to use the '-M' option name and not '--smpiarg' since
  # 'trilinos_jsrun' has special logic
fi

# NOTE: We used to check for the launch node but at one point that changed
# from 'vortex59' to 'vortex5' without warning.  That caused all of the tests
# run with 'trilinos_jsrun' to fail on 2020-08-11 so we got no test results.
# Therefore, we will not be checking for running on the launch node anymore in
# order to avoid having all of the testing break when they change the launch
# node again.  Therefore, we also removed checks for the login node as well.
# This makes these scripts more robust to changes on 'vortex'.


#
# Define functions for running on compute nodes
#


function atdm_ats2_get_allocated_compute_node_name() {
  if [[ "${LSB_HOSTS}" != "" ]] ; then
    echo ${LSB_HOSTS} | cut -d' ' -f2
  else
     echo
  fi
}
export -f atdm_ats2_get_allocated_compute_node_name


#
# All done!
#

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
