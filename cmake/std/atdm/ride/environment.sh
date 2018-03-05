################################################################################
#
# Set up env on ride/white for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
################################################################################

module purge

export OMPI_CXX=
export BUILD_COUNT=128

echo "Using white/ride compiler stack $ATDM_CONFIG_COMPILER to build $ATDM_CONFIG_BUILD_TYPE code with Kokkos node type $ATDM_CONFIG_NODE_TYPE"

module load ninja/1.7.2

export ATDM_CONFIG_KOKKOS_ARCH=Power8
if [ "$ATDM_CONFIG_COMPILER" == "GNU" ]; then
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    export OMPI_CXX=`which g++`
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp"
elif [ "$ATDM_CONFIG_COMPILER" == "INTEL" ]; then
    echo "Intel compiler not supported"
    return
elif [ "$ATDM_CONFIG_COMPILER" == "CUDA" ]; then
    export ATDM_CONFIG_KOKKOS_ARCH=Kepler37
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    export ATDM_CONFIG_LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export ATDM_CONFIG_BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp"
    export OMPI_CXX=$ATDM_CONFIG_TRILNOS_DIR/packages/kokkos/config/nvcc_wrapper
    if [ ! -x "$OMPI_CXX" ]; then
	export OMPI_CXX=`which nvcc_wrapper`
    fi
    if [ ! -x "$OMPI_CXX" ]; then
        echo "No nvcc_wrapper found"
        exit 1
    fi
    export OMPI_CC=`which gcc`
    export OMPI_FC=`which gfortran`
    export ATDM_CONFIG_USE_CUDA=ON
    export CUDA_LAUNCH_BLOCKING=1
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
else
    echo "No valid compiler found"

fi

export ATDM_CONFIG_MPI_POST_FLAG="-map-by;socket:PE=8;--oversubscribe"

export ATDM_CONFIG_USE_HWLOC=OFF

export ATDM_CONFIG_HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export ATDM_CONFIG_NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

module swap yamlcpp/0.3.0 yaml-cpp/20170104 
if [ $? ]; then module load  yaml-cpp/20170104; fi

# If you are running OpenMP
export OMP_NUM_THREADS=2

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`
