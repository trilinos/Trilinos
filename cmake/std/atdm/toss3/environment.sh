################################################################################
#
# Set up env on toss3 (chama/serrano) for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

echo "Using chama compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

# there does not appear to be a ninja module ontoss3 so turn off ninja
export ATDM_CONFIG_USE_MAKEFILES=ON
export ATDM_CONFIG_BUILD_COUNT=16
export ATDM_CONFIG_USE_NINJA=OFF

module purge
. /projects/sems/modulefiles/utils/sems-modules-init.sh
module load sems-env
module load cmake/3.5


if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

#export Kokkos_Arch=BDW
if [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load sems-python/2.7.9
    module load sems-intel/17.0.0
    module load sems-openmpi/1.10.5
    module load sems-hdf5/1.8.12/parallel
    module load sems-netcdf/4.4.1/exo_parallel 
    module load sems-yaml_cpp/0.5.3/base 
    module load sems-boost/1.59.0/base  
    module load intel/17.0.4.196
    module load mkl/17.0.4.196
    export BOOST_ROOT=$SEMS_BOOST_ROOT
    export HDF5_ROOT=$SEMS_HDF5_ROOT
    export NETCDF_ROOT=$SEMS_NETCDF_ROOT
    export YAMLCPP_ROOT=$SEMS_YAML_CPP_ROOT
    export OMPI_CXX=`which icpc`
    export OMPI_CC=`which icc`
    export OMPI_FC=`which ifort`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
else
    echo "No valid compiler found"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${SEMS_HDF5_ROOT}/lib;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${SEMS_BOOST_ROOT}/lib;-L${SEMS_NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${SEMS_HDF5_ROOT}/lib;${SEMS_BOOST_ROOT}/lib/libboost_program_options.a;${SEMS_BOOST_ROOT}/lib/libboost_system.a;${SEMS_NETCDF_ROOT}/lib/libnetcdf.a;${SEMS_NETCDF_ROOT}/lib/libpnetcdf.a;${SEMS_HDF5_ROOT}/lib/libhdf5_hl.a;${SEMS_HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lcurl"

# not sure what below does.  It was in the original environment script
#unset ATTB_ENV

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=16;--oversubscribe"

# Set the default compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif77
export F90=mpif90

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
