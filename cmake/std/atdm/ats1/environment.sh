################################################################################
#
# Set up env on ats1 for ATDM builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL-19.0.4_MPICH-7.7.6
fi

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

# Common sparc tpl path values
sparc_tpl_prefix_path="/usr/projects/sparc/tpls"
system_name="ats1"
node_arch=""

if [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "HSW" ]]; then
  node_arch="hsw"
  module unload craype-mic-knl
  module load craype-haswell
  # HSW nodes have 64 virtual cores.
  # Allow no more than 8 virtual cores per task.
  export ATDM_CONFIG_MPI_POST_FLAGS="-c 4"
  # If we have 1 MPI rank per srun command and 1 cpu per task, we can run up to 32 2-threaded tests "in parallel" on virtual cores.
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=1
elif [[ "$ATDM_CONFIG_KOKKOS_ARCH" == "KNL" ]]; then
  node_arch="knl"
  module unload craype-haswell
  module load craype-mic-knl
  export SLURM_TASKS_PER_NODE=16
  export OMP_PLACES=threads
  export OMP_PROC_BIND=spread
  # Ensure that no more than 8 tasks, per srun command, are launched.
  export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node=8"
  # KNL nodes have 272 virtual cores.
  # Allow now more than 34 virtual cores per task.
  export ATDM_CONFIG_MPI_POST_FLAGS="--hint=nomultithread;-c 4"
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=1
  export ATDM_CONFIG_SBATCH_EXTRA_ARGS="-p knl -C cache --hint=multithread"
  export ATDM_CONFIG_BUILD_COUNT=8
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not"
  echo "*** a supported combination on this system!"
  echo "*** Combinations that are supported: "
  echo "*** > intel-19.0.4_mpich-7.7.6"
  echo "*** > intel-18.0.5_mpich-7.7.6"
  echo "***"
  return
fi

echo "Using ats1 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH"

# Exclude bad nodes.
export ATDM_CONFIG_SBATCH_EXTRA_ARGS="$ATDM_CONFIG_SBATCH_EXTRA_ARGS --exclude=nid00021,nid00020"

export OMP_NUM_THREADS=2

# Common modules and paths
module load cmake/3.14.6
module load git/2.19.1
export PATH=/projects/netpub/atdm/ninja-1.8.2/bin:$PATH
module unload cray-mpich
module load cray-mpich/7.7.6
module unload cray-libsci
module unload gcc

if [[ "$ATDM_CONFIG_COMPILER" == "INTEL-19.0.4"* ]]; then
  module load gcc/8.2.0
  module unload intel
  module load intel/19.0.4
  #sparc_tpl_ext=${system_name}-${node_arch}_intel-19.0.4
  #sparc_tpl_mpi_ext=${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
elif [[ "$ATDM_CONFIG_COMPILER" == "INTEL-18.0.5"* ]]; then
  module unload intel
  module load intel/18.0.5
  #sparc_tpl_ext=${system_name}-${node_arch}_intel-18.0.5
  #sparc_tpl_mpi_ext=${system_name}-${node_arch}_intel-18.0.5_mpich-7.7.6
else
  echo
  echo "***"
  echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER and KOKKOS_ARCH=$ATDM_CONFIG_KOKKOS_ARCH is not"
  echo "*** a supported combination on this system!"
  echo "*** Combinations that are supported: "
  echo "*** > intel-19.0.4_mpich-7.7.6_openmp-HSW"
  echo "*** > intel-19.0.4_mpich-7.7.6_openmp-KNL"
  echo "*** > intel-18.0.5_mpich-7.7.6_openmp-HSW"
  echo "*** > intel-18.0.5_mpich-7.7.6_openmp-KNL"
  echo "***"
  return
fi

export PATH=/usr/projects/hpcsoft/cle6.0/common/intel-clusterstudio/2019.4-068/compilers_and_libraries_2019/linux/mkl/bin:$PATH
export LD_LIBRARY_PATH=/usr/projects/hpcsoft/cle6.0/common/intel-clusterstudio/2019.4-068/compilers_and_libraries_2019/linux/mkl/lib/intel64:$LD_LIBRARY_PATH
# This linker path is needed for cmake's C, CXX, and Ftn compiler tests. This linker path is not shown
# anywhere in SPARC's environment but is needed in ATDM's environment. Cmake's compiler tests are
# ignoring LD_{CRAY,CRAYPAT,}_LIBRARY_PATH for some reason, so we add this path via LDFLAGS via
# s.t. it is added to CMAKE_EXE_LINKER_FLAGS before trilinos's cmake config invokes PROJECT().
export LDFLAGS="-L/opt/gcc/8.2.0/snos/lib/gcc/x86_64-suse-linux/8.2.0/ $LDFLAGS"

export CBLAS_ROOT=/usr/projects/hpcsoft/cle6.0/common/intel-clusterstudio/2019.4-068/compilers_and_libraries_2019/linux
export COMPILER_ROOT=/usr/projects/hpcsoft/cle6.0/common/intel-clusterstudio/2019.4-068/compilers_and_libraries_2019/linux
export MPI_ROOT=/opt/cray/pe/mpt/7.7.6/gni/mpich-intel/16.0

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

# Boost 1.65.1 settings
export BOOST_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/boost-1.72.0/00000000/${system_name}-${node_arch}_intel-19.0.4

# Hdf5 1.10.5 settings
export HDF5_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/hdf5-1.10.5/00000000/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
#export ATDM_CONFIG_HDF5_NO_SYSTEM_PATHS=ON

# BINUTILS settings
export BINUTILS_ROOT="/usr"
export ATDM_CONFIG_BINUTILS_LIBS="${BINUTILS_ROOT}/lib64/libbfd.a;-lz;${BINUTILS_ROOT}/lib64/libiberty.a"

# Cgns settings
export CGNS_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/cgns-c09a5cd/d313cc2f822078e47c7dbdee074ecb0431e573eb/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
export ATDM_CONFIG_CGNS_LIBRARY_NAMES="cgns"

# Pnetcdf 1.10.0 settings
export PNETCDF_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/pnetcdf-1.12.1/6144dc67b2041e4093063a04e89fc1e33398bd09/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6

# Netcdf 4.7.0 settings
export NETCDF_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/netcdf-4.7.0/24baa07a3fa1ff9dbc8e70dc591ebbdec56783b2/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib64;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lm"

# Libhio 1.4.1.2 settings
export LIBHIO_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/libhio-1.4.1.2/00000000/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6

# Metis 4.0.3 settings
export METIS_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/parmetis-4.0.3/00000000/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6

# Parmetis 4.0.3 settings
export PARMETIS_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/parmetis-4.0.3/00000000/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
# Force linker order!
export ATDM_CONFIG_PARMETIS_LIBS="${PARMETIS_ROOT}/lib/libparmetis.a;${METIS_ROOT}/lib/libmetis.a"

# Superludist 5.4.0 settings
export SUPERLUDIST_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/superlu_dist-5.4.0/a3121eaff44f7bf7d44e625c3b3d2a9911e58876/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS="${SUPERLUDIST_ROOT}/include"
export ATDM_CONFIG_SUPERLUDIST_LIBS="${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a"

# Sgm 20.23 settings
export SGM_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/sgm-20.23/00000000/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6

# Euclid 20.23 settings
export EUCLID_ROOT=${sparc_tpl_prefix_path}/${system_name}-${node_arch}/euclid-20.23/8b68b12f72b59648c9a0a962a6d55ea978199860/${system_name}-${node_arch}_intel-19.0.4_mpich-7.7.6

# Define and export atdm_run_script_on_compute_node
source $ATDM_SCRIPT_DIR/common/define_run_on_slurm_compute_node_func.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
