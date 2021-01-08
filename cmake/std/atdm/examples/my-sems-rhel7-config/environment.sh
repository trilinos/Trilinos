# Example of how to set up your own custom ATDM Trilinos configuration on some
# machine.  This example only works on SEMS RHE6 systems.

echo "Using my custom SEMS RHEL7 env for compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

if [[ "${ATDM_CONFIG_COMPILER}" != "GNU" ]] ; then
  echo "ATDM_CONFIG_COMPILER=${ATDM_CONFIG_COMPILER} not supported!  Only 'gnu' supported currently!"
fi

# 1) Load modules, set paths, etc.

module purge
module load sems-env

module load sems-git/2.10.1
module load sems-cmake/3.17.1
module load sems-ninja_fortran/1.10.0

module load sems-gcc/7.2.0
module load sems-openmpi/1.10.1
module load sems-netcdf/4.7.3/parallel
module load sems-hdf5/1.10.6/parallel
module load sems-zlib/1.2.8/base
module load sems-boost/1.59.0/base
module unload sems-python/2.7.9
module load sems-superlu/4.3/base

# 2) Set up env vars that will get picked up by ATDMDevEnvSettings.cmake when
# configuring Trilinos.

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=20
export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=1
unset ATDM_CONFIG_KOKKOS_ARCH # Or set to the arch you know you have!

# MPI Stuff
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`
export OMPI_CXX=`which g++`
export OMPI_CC=`which gcc`
export OMPI_FC=`which gfortran`
export ATDM_CONFIG_MPI_PRE_FLAGS="--bind-to;none"

# OpenMP stuff
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=false
unset OMP_PLACES

# TPL stuff
export ATDM_CONFIG_USE_HWLOC=OFF
export LAPACK_ROOT=/usr/lib64
export ATDM_CONFIG_LAPACK_LIBS="-L${LAPACK_ROOT};-llapack"
export ATDM_CONFIG_BLAS_LIBS="-L${LAPACK_ROOT}/lib;-lblas"
export ZLIB_ROOT=${SEMS_ZLIB_ROOT}
export BOOST_ROOT=${SEMS_BOOST_ROOT}
export HDF5_ROOT=${SEMS_HDF5_ROOT}
export NETCDF_ROOT=${SEMS_NETCDF_ROOT}

# Done
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
