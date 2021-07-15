################################################################################
#
# Set up env on blake (intel skylake) for ATMD builds of Trilinos
#
# This source script gets the settings from the ATDM_CONFIG_BUILD_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
    export ATDM_CONFIG_COMPILER=ONEAPI_2021.2.0
elif [ "$ATDM_CONFIG_COMPILER" == "ONEAPI-2021.2.0" ] ; then
    export ATDM_CONFIG_COMPILER=ONEAPI_2021.2.0
elif [ "$ATDM_CONFIG_COMPILER" == "ONEAPI-2021.1.1" ] ; then
    export ATDM_CONFIG_COMPILER=ONEAPI_2021.1.1
fi

echo "Using blake compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=12
# export ATDM_CONFIG_CMAKE_JOB_POOL_LINK=2
# NOTE: Above, currently setting CMAKE_JOB_POOL_LINK results in build
# failures with Ninja.  See https://gitlab.kitware.com/snl/project-1/issues/60

module purge

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
    export OMP_NUM_THREADS=2
else
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
fi

if [ "$ATDM_CONFIG_COMPILER" == "ONEAPI" ]; then
    module load devpack/20210420/openmpi/4.0.5/intel/oneapi/2021.2.0

    export OMPI_CXX=`which icpx`
    export OMPI_CC=`which icx`
    export OMPI_FC=`which ifx`
    export ATDM_CONFIG_LAPACK_LIBS="-mkl"
    export ATDM_CONFIG_BLAS_LIBS="-mkl"
elif [ "$ATDM_CONFIG_COMPILER" == "ONEAPI_2021.2.0" ]; then
    module load devpack/20210420/openmpi/4.0.5/intel/oneapi/2021.2.0

    export OMPI_CXX=`which icpx`
    export OMPI_CC=`which icx`
    export OMPI_FC=`which ifx`
    export ATDM_CONFIG_LAPACK_LIBS="-mkl"
    export ATDM_CONFIG_BLAS_LIBS="-mkl"
elif [ "$ATDM_CONFIG_COMPILER" == "ONEAPI_2021.1.1" ]; then
    module load devpack/20210310/openmpi/4.0.5/intel/oneapi/2021.1.1

    export OMPI_CXX=`which icpx`
    export OMPI_CC=`which icx`
    export OMPI_FC=`which ifx`
    export ATDM_CONFIG_LAPACK_LIBS="-mkl"
    export ATDM_CONFIG_BLAS_LIBS="-mkl"
else
    echo
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${NETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${NETCDF_ROOT}/lib64/libnetcdf.a;${PNETCDF_ROOT}/lib64/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS};-lcurl"


# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_EXEC=mpirun
export ATDM_CONFIG_MPI_PRE_FLAGS=""
export ATDM_CONFIG_MPI_EXEC_NUMPROCS_FLAG=--np

# Set the default compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif77
export F90=mpif90

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
