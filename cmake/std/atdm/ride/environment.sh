#!/bin/bash
module purge

export OMPI_CXX=
BUILD_COUNT=128

if [ -z "$INSTALL" ]; then export INSTALL=OFF; fi
if [ -z "$PROJECT" ]; then PROJECT=PANZER; fi
if [ -z "$TEST" ]; then export TEST=OFF; fi

export BUILD_DIR=build_${PROJECT}_${COMPILER}_${BUILD_TYPE}_${NODE_TYPE}
export INSTALL_DIR=/home/projects/empire/ride/install-for-${PROJECT}/$COMPILER-$BUILD_TYPE-$NODE_TYPE
export Trilinos_INSTALL_DIR=`echo $INSTALL_DIR | sed s/${PROJECT}/PANZER/`
export PATH=$Trilinos_INSTALL_DIR/bin:$PATH

if [ $INSTALL == "ON" ]; then
    if [[ ! -e $INSTALL_DIR ]]; then mkdir -p $INSTALL_DIR; fi
    if [[ ! -e $INSTALL_DIR ]]; then mkdir -p $INSTALL_DIR; fi
    echo "#!/bin/bash" >$INSTALL_DIR/configure.sh
    echo "JOB_NAME="$JOB_NAME >>$INSTALL_DIR/configure.sh
    echo export Trilinos_INSTALL_DIR=$INSTALL_DIR>>$INSTALL_DIR/configure.sh
    echo . $INSTALL_DIR/set_build_options.sh>>$INSTALL_DIR/configure.sh
fi


export Kokkos_Arch=Power8
if [ "$COMPILER" == "GNU" ]; then
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    export OMPI_CXX=`which g++`
    export LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp"
elif [ "$COMPILER" == "INTEL" ]; then
    echo "Intel compiler not supported"
    exit 1

elif [ "$COMPILER" == "CUDA" ]; then
    export Kokkos_Arch=Kepler37
    module load devpack/openmpi/1.10.4/gcc/5.4.0/cuda/8.0.44
    export LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran;-lgomp"
    export BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran;-lgomp"
    export OMPI_CXX=$WORKSPACE/Trilinos/packages/kokkos/config/nvcc_wrapper
    if [ $INSTALL == "ON" ]; then
        if [[ ! -e $INSTALL_DIR/bin ]]; then mkdir -p $INSTALL_DIR/bin; fi
        cp $WORKSPACE/Trilinos/packages/kokkos/config/nvcc_wrapper $INSTALL_DIR/bin
    fi
    if [ ! -x "$OMPI_CXX" ]; then
	export OMPI_CXX=`which nvcc_wrapper`
        if [ ! -x "$OMPI_CXX" ]; then
            export OMPI_CXX=$INSTALL_DIR/bin/nvcc_wrapper 
        fi
    fi
    if [ ! -x "$OMPI_CXX" ]; then
        echo "No nvcc_wrapper found"
        exit 1
    fi
    export USE_CUDA=ON
    export CUDA_LAUNCH_BLOCKING=1
    export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
else
    echo "No valid compiler found"

fi
export MPI_POST_FLAG="-map-by;socket:PE=8;--oversubscribe"

export USE_HWLOC=OFF

export HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"

module swap yamlcpp/0.3.0 yaml-cpp/20170104 
if [ $? ]; then module load  yaml-cpp/20170104; fi
export COMPILER
unset ATTB_ENV
# Set the default compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif77
export F90=mpif90
