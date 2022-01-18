################################################################################
#
# Set up env on toss3 (chama and serrano) for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

echo "Using $ATDM_CONFIG_SYSTEM_NAME toss3 compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_ENABLE_SPARC_SETTINGS=ON
export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=8
# export ATDM_CONFIG_CMAKE_JOB_POOL_LINK=2
# NOTE: Above, currently setting CMAKE_JOB_POOL_LINK results in a build
# failures with Ninja.  See https://gitlab.kitware.com/snl/project-1/issues/60

if [[ "${ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE}" != "1" ]] ; then
  # We do this twice since sems modules are wacked and we get errors to the
  # screen on a purge The second purge will catch any real errors with purging
  # ...
  module purge &> /dev/null
  module purge
else
  echo "NOTE: ATDM_CONFIG_DONT_LOAD_SPARC_MODULES_PLEASE=1 is set so using pre-loaded sparc-dev module!"
fi

. /projects/sems/modulefiles/utils/sems-archive-modules-init.sh
module load sems-archive-env
module load sems-archive-git/2.10.1

# Common paths and modules for both intel-1{8,9}
module load sems-archive-cmake/3.19.1

atdm_config_load_sparc_dev_module sparc-dev/intel-19.0.4_openmpi-4.0.3

if [ "$ATDM_CONFIG_COMPILER" == "INTEL-19.0.4_OPENMPI-4.0.3" ]; then
  # Correct module already loaded above
  :
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL-18.0.2_OPENMPI-4.0.3" ]; then
  module swap intel/19.0 intel/18.0.2.199 &> /dev/null
  atdm_remove_dirs_from_path "/projects/sparc/tools/vvt"
else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

module load sems-archive-ninja_fortran/1.8.2

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
  unset OMP_PLACES
  unset OMP_PROC_BIND
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

export OMPI_CXX=`which icpc`
export OMPI_CC=`which icc`
export OMPI_FC=`which ifort`
export ATDM_CONFIG_LAPACK_LIBS="-mkl"
export ATDM_CONFIG_BLAS_LIBS="-mkl"

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"
export ATDM_CONFIG_SUPERLUDIST_INCLUDE_DIRS=${SUPERLUDIST_ROOT}/include
export ATDM_CONFIG_SUPERLUDIST_LIBS=${SUPERLUDIST_ROOT}/lib64/libsuperlu_dist.a
export ATDM_CONFIG_BINUTILS_LIBS="/usr/lib64/libbfd.so;/usr/lib64/libiberty.a"

# not sure what below does.  It was in the original environment script
#unset ATTB_ENV

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_EXEC=srun
export ATDM_CONFIG_MPI_PRE_FLAGS="--mpi=pmi2;--ntasks-per-node;36"
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=--ntasks

# Define function atdm_run_script_on_compute_node
source $ATDM_SCRIPT_DIR/common/define_run_on_slurm_compute_node_func.sh

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
