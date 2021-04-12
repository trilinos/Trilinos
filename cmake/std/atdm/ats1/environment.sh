################################################################################
#
# Set up env on ats1 for ATDM builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################
if [ "$ATDM_CONFIG_KOKKOS_ARCH" == "DEFAULT" ] ; then
  export ATDM_CONFIG_KOKKOS_ARCH=HSW
fi

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON

# Use srun to run mpi jobs
export ATDM_CONFIG_MPI_EXEC="/opt/slurm/bin/srun"

# srun does not accept "-np" for # of processors
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="--ntasks"
# Use pmi2 as MPI type.
# Ensure that no more than 36 tasks, per srun command, are launched.
export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node;36"

export ATDM_CONFIG_SBATCH_DEFAULT_ACCOUNT=IGNORED

# Assume we are building on a haswell compute node with 64 virtual cores
export ATDM_CONFIG_BUILD_COUNT=32
# NOTE: Above, we were getting out-of-memory errors when trying to build with
# 64 processes so we reduced this to just 32 to try to avoid these.  (See
# ATDV-361)

if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW" ]]; then
  atdm_config_load_sparc_dev_module sparc-dev/intel-19.0.4_mpich-7.7.15_hsw
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  export OMP_NUM_THREADS=2
  unset OMP_PLACES
elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "KNL" ]]; then
  atdm_config_load_sparc_dev_module sparc-dev/intel-19.0.4_mpich-7.7.15_knl
  export SLURM_TASKS_PER_NODE=16
  # Ensure that no more than 8 tasks, per srun command, are launched.
  export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node=8"
  # KNL nodes have 272 virtual cores.
  # Allow now more than 34 virtual cores per task.
  export ATDM_CONFIG_MPI_POST_FLAGS="--hint=nomultithread;-c 8"
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=1
  export ATDM_CONFIG_SBATCH_EXTRA_ARGS="-p knl -C cache --hint=multithread"
  export ATDM_CONFIG_BUILD_COUNT=8
  export OMP_NUM_THREADS=8
  unset OMP_PLACES
  unset OMP_PROC_BIND
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not"
  echo "*** a supported combination on this system!"
  echo "*** Combinations that are supported: "
  echo "*** > intel-19.0.4_mpich-7.7.15"
  echo "***"
  return
fi

if [[ "${SPARC_MODULE}" == "" ]] ; then
  echo "No SPARC_MODULE loaded so exiting!"
  return
fi

echo "Using ats1 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

# Exclude bad nodes.
export ATDM_CONFIG_SBATCH_EXTRA_ARGS="$ATDM_CONFIG_SBATCH_EXTRA_ARGS --exclude=nid00021,nid00020"

# Common modules and paths
export PATH=/projects/netpub/atdm/ninja-1.8.2/bin:$PATH

# This linker path is needed for cmake's C, CXX, and Ftn compiler tests. This linker path is not shown
# anywhere in SPARC's environment but is needed in ATDM's environment. Cmake's compiler tests are
# ignoring LD_{CRAY,CRAYPAT,}_LIBRARY_PATH for some reason, so we add this path via LDFLAGS via
# s.t. it is added to CMAKE_EXE_LINKER_FLAGS before trilinos's cmake config invokes PROJECT().
export LDFLAGS="-L/opt/gcc/8.3.0/snos/lib/gcc/x86_64-suse-linux/8.3.0/ -lpthread $LDFLAGS"

# 2021-01-09 -- Resolve linker errors in Seacas and Percept
export LDFLAGS="-L${MPI_ROOT}/lib -lmpich -lrt ${ATP_INSTALL_DIR}/lib/libAtpSigHandler.a ${ATP_INSTALL_DIR}/lib/libbreakpad_client_nostdlib.a $LDFLAGS"

export ATDM_CONFIG_TPL_FIND_SHARED_LIBS=OFF
export ATDM_CONFIG_Trilinos_LINK_SEARCH_START_STATIC=ON
export ATDM_CONFIG_Trilinos_ENABLE_SECONDARY_TESTED_CODE=OFF
export ATDM_CONFIG_CMAKE_SKIP_INSTALL_RPATH=ON

# Set MPI wrappers
export MPICXX=`which CC`
export MPICC=`which cc`
export MPIF90=`which ftn`

# Set common default compilers
export CC=${MPICC}
export CXX=${MPICXX}
export F77=${MPIF90}
export FC=${MPIF90}
export F90=${MPIF90}

# Anasazi settings
export ATDM_CONFIG_Anasazi_ENABLE_RBGen=OFF

# Lapack (intel) settings
export ATDM_CONFIG_LAPACK_LIBS="-L${CBLAS_ROOT}/mkl/lib/intel64;-L${CBLAS_ROOT}/compiler/lib/intel64;-mkl;-lmkl_intel_lp64;-lmkl_intel_thread;-lmkl_core;-liomp5"

# Blas (intel) settings
export ATDM_CONFIG_BLAS_LIBS="-L${CBLAS_ROOT}/mkl/lib/intel64;-L${CBLAS_ROOT}/compiler/lib/intel64;-mkl;-lmkl_intel_lp64;-lmkl_intel_thread;-lmkl_core;-liomp5"

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
#export ATDM_CONFIG_HDF5_NO_SYSTEM_PATHS=ON

# BINUTILS settings
export BINUTILS_ROOT="/usr"
export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib64/libbfd.a;-lz;${BINUTILS_ROOT}/lib64/libiberty.a"

export ATDM_CONFIG_CGNS_LIBRARY_NAMES="cgns"

export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib64;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lm"

# Force linker order!
export ATDM_CONFIG_PARMETIS_LIBS="${PARMETIS_ROOT}/lib/libparmetis.a;${METIS_ROOT}/lib/libmetis.a"

export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS="${SUPERLUDIST_ROOT}/include"
export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a"

# Define and export atdm_run_script_on_compute_node
source $ATDM_SCRIPT_DIR/common/define_run_on_slurm_compute_node_func.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
