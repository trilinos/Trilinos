################################################################################
#
# Set up env on ride/white for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

if [ "$ATDM_CONFIG_COMPILER" == "DEFAULT" ] ; then
  export ATDM_CONFIG_COMPILER=GNU
fi

echo "Using white/ride compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

export ATDM_CONFIG_USE_NINJA=ON
export ATDM_CONFIG_BUILD_COUNT=128
# NOTE: Above settings are used for running on a single rhel7F (Firestone,
# Dual-Socket POWER8, 8 cores per socket, K80 GPUs) node.

module purge

if [ "$ATDM_CONFIG_NODE_TYPE" == "OPENMP" ] ; then
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
  export OMP_NUM_THREADS=2
else
  export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=32
fi

export ATDM_CONFIG_KOKKOS_ARCH=Power8,Kepler37
if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    module swap openblas/0.2.19/gcc/5.4.0 netlib/3.8.0/gcc/5.4.0
    export OMPI_CXX=`which g++`
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"
elif [ "$ATDM_CONFIG_COMPILER" == "CUDA" ]; then
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=16
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    module swap openblas/0.2.19/gcc/5.4.0 netlib/3.8.0/gcc/5.4.0
    export OMPI_CXX=$ATDM_CONFIG_TRILNOS_DIR/packages/kokkos/bin/nvcc_wrapper
    if [ ! -x "$OMPI_CXX" ]; then
        echo "No nvcc_wrapper found"
        return
    fi
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp;-lm"
    export ATDM_CONFIG_USE_CUDA=ON
    export CUDA_LAUNCH_BLOCKING=1
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
    export ATDM_CONFIG_CTEST_PARALLEL_LEVEL=8
    # Avoids timeouts due to not running on separate GPUs (see #2446)
else
    echo "***"
    echo "*** ERROR: COMPILER=$ATDM_CONFIG_COMPILER is not supported on this system!"
    echo "***"
    return
fi

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

module swap yamlcpp/0.3.0 yaml-cpp/20170104 
if [ $? ]; then module load  yaml-cpp/20170104; fi

module load binutils/2.27.0
# NOTE: Above we need to use updated 'ar' from binutils 2.27 to handle big
# object files with 'ar' that otherwise can cause "File truncated" failures (see #3069)

# Use manually installed cmake and ninja to try to avoid module loading
# problems (see TRIL-208)
export PATH=/ascldap/users/rabartl/install/white-ride/cmake-3.11.2/bin:/ascldap/users/rabartl/install/white-ride/ninja-1.8.2/bin:$PATH

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export ATDM_CONFIG_MPI_PRE_FLAGS="--mca;orte_abort_on_non_zero_status;0"
export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=4"

export ATDM_CONFIG_COMPLETED_ENV_SETUP=TRUE
