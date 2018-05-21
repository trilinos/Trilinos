################################################################################
#
# Set up env on mutrino for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=INTEL
fi

echo "Using mutrino compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=OFF
export ATDM_CONFIG_BUILD_COUNT=8
export OMP_NUM_THREADS=2

# Don't purge as this removes the Cray default environment
#module purge
module use /projects/EMPIRE/mutrino/tpls/hsw/modulefiles

# Use srun to run mpi jobs
export ATDM_CONFIG_MPI_EXEC="/opt/slurm/bin/srun"

# srun does not accept "-np" for # of processors
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG="-n"
export ATDM_CONFIG_KOKKOS_ARCH=HSW
export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16

if [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load devpack/20180124/cray/7.6.2/intel/17.0.4
    module load gcc/4.9.3
    module load cmake/3.9.0
    export MPICXX=`which CC`
    export MPICC=`which cc`
    export MPIF90=`which ftn`

#    # Cray provides differently named wrappers
#    export CXX=`which CC`
#    export CC=`which cc`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
else
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl;-lm"

export ATDM_CONFIG_MPI_POST_FLAG="-c 4"
export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
