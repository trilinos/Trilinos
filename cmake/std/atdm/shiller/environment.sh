################################################################################
#
# Set up env on shiller/hansen for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

module purge

export OMPI_CXX=

export ATDM_CONFIG_SYSTEM_CDASH_SITE=hansen/shiller
export ATDM_CONFIG_USE_MAKEFILES=ON
export ATDM_CONFIG_BUILD_COUNT=32
export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
export OMP_NUM_THREADS=2

echo "Using hansen/shiller compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

module load ninja/1.7.2

export ATDM_CONFIG_KOKKOS_ARCH=HSW
if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=`which g++`
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    module load devpack/openmpi/2.1.1/intel/17.4.196/cuda/none
    export OMPI_CXX=`which icpc`
    export OMPI_CC=`which icc`
    export OMPI_FC=`which ifort`
    export ATDM_CONFIG_LAPACK_LIB="-mkl"
    export ATDM_CONFIG_BLAS_LIB="-mkl"
elif [ "$ATDM_CONFIG_COMPILER" == "CUDA" ]; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=$ATDM_CONFIG_TRILNOS_DIR/packages/kokkos/config/nvcc_wrapper 
    if [ ! -x "$OMPI_CXX" ]; then
	export OMPI_CXX=`which nvcc_wrapper`
    fi
    if [ ! -x "$OMPI_CXX" ]; then
        echo "No nvcc_wrapper found"
        return  # Not 'exit' since this is sourced, not called!
    fi
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_USE_CUDA=ON
    export CUDA_LAUNCH_BLOCKING=1
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
else
    echo "No valid compiler found"
fi

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_USE_HWLOC=OFF
export ATDM_CONFIG_HWLOC_LIBS=-lhwloc

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;-lhdf5_hl;-lhdf5;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${ATDM_CONFIG_HDF5_LIBS}"

export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=16;--oversubscribe"

# Just in case the user runs with an OpenMP build
export OMP_NUM_THREADS=2
