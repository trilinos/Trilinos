################################################################################
#
# Set up env on shiller/hansen for ATMD builds of Trilinos
#
# This source script gets the settings from the JOB_NAME var.
#
#  
#
################################################################################

# Assert this script is sourced, not run!
called=$_
if [ "$called" == "$0" ] ; then
  echo "This script '$0' is being called.  Instead, it must be sourced!"
  exit 1
fi

# Assert that JOB_NAME is set!
if [ -z "$JOB_NAME" ] ; then
  echo "Error, must set JOB_NAME in env!"
  return
fi

# Return the absoute directory of some relative directory path.
#
# This uses a temp shell to cd into the directory and then uses pwd to get the
# path.
function get_abs_dir_path() {
  [ -z "$1" ] && { pwd; return; }
  (cd -P -- "$1" && pwd)
}

# Get the base dir for the sourced script to find the base of Trilinos
_SCRIPT_DIR=`echo $BASH_SOURCE | sed "s/\(.*\)\/.*\.sh/\1/g"`
#echo "_SCRIPT_DIR = '$_SCRIPT_DIR'"
TRILNOS_DIR=`get_abs_dir_path $_SCRIPT_DIR/../../../..`
echo "TRILNOS_DIR = '$TRILNOS_DIR'"

# Purge all of the existing modules first
module purge

export BUILD_COUNT=32
export OMPI_CXX=

# Set the defaults
export COMPILER=GNU
export BUILD_TYPE=DEBUG
export USE_OPENMP=OFF
export USE_CUDA=OFF
export USE_PTHREADS=OFF

if [[ $JOB_NAME == *"Trilinos"* ]]; then export PROJECT=PANZER; fi
if [[ $JOB_NAME == *"PIC"* ]]; then export PROJECT=EMPIRE-PIC; fi
if [[ $JOB_NAME == *"openmp"* ]]; then export USE_OPENMP=ON; fi
if [[ $JOB_NAME == *"opt"* ]]; then export BUILD_TYPE=RELEASE; fi
if [[ $JOB_NAME == *"gnu"* ]]; then COMPILER=GNU; fi
if [[ $JOB_NAME == *"intel"* ]]; then COMPILER=INTEL; fi
if [[ $JOB_NAME == *"cuda"* ]]; then 
    export USE_CUDA=ON; 
    COMPILER=CUDA; 
fi

if [ -z "$INSTALL" ]; then export INSTALL=OFF; fi
if [ -z "$PROJECT" ]; then PROJECT=PANZER; fi
if [ -z "$TEST" ]; then export TEST=OFF; fi

if [ "$USE_OPENMP" == "ON" ]; then 
    NODE_TYPE=OPENMP
elif [ "$USE_PTHREADS" == "ON" ]; then 
    NODE_TYPE=THREAD
elif [ "$USE_CUDA" == "ON" ]; then 
    export COMPILER=CUDA
    NODE_TYPE=CUDA
else
    NODE_TYPE=SERIAL
fi
export NODE_TYPE

export BUILD_DIR=build_${PROJECT}_${COMPILER}_${BUILD_TYPE}_${NODE_TYPE}
export INSTALL_DIR=/home/projects/empire/install-for-${PROJECT}/$COMPILER-$BUILD_TYPE-$NODE_TYPE
export Trilinos_INSTALL_DIR=`echo $INSTALL_DIR | sed s/${PROJECT}/PANZER/`
export PATH=$Trilinos_INSTALL_DIR/bin:$PATH

if [ $INSTALL == "ON" ]; then
    if [[ ! -e $INSTALL_DIR ]]; then mkdir -p $INSTALL_DIR; fi
    echo "#!/bin/bash" >$INSTALL_DIR/configure.sh
    echo "JOB_NAME="$JOB_NAME >>$INSTALL_DIR/configure.sh
    echo export Trilinos_INSTALL_DIR=$INSTALL_DIR>>$INSTALL_DIR/configure.sh
fi

echo using compiler stack $COMPILER to build $BUILD_TYPE code

if [ "$COMPILER" == "GNU" ]; then
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=`which g++`
    export LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
elif [ "$COMPILER" == "INTEL" ]; then
    module load devpack/openmpi/2.1.1/intel/17.4.196/cuda/none
    export OMPI_CXX=`which icpc`
    export OMPI_CC=`which icc`
    export OMPI_FC=`which ifort`
    export LAPACK_LIB="-mkl"
    export BLAS_LIB="-mkl"
elif [ "$COMPILER" == "CUDA" ]; then
    module load devpack/openmpi/2.1.1/gcc/4.9.3/cuda/8.0.61
    export OMPI_CXX=$TRILNOS_DIR/packages/kokkos/config/nvcc_wrapper 
    if [ $INSTALL == "ON" ]; then
        if [[ ! -e $INSTALL_DIR/bin ]]; then mkdir -p $INSTALL_DIR/bin; fi
        cp $TRILNOS_DIR/packages/kokkos/config/nvcc_wrapper $INSTALL_DIR/bin
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
    export LAPACK_LIB="-L${LAPACK_ROOT}/lib;-llapack;-lgfortran"
    export BLAS_LIB="-L${BLAS_ROOT}/lib;-lblas;-lgfortran"
else
    echo "No valid compiler found"
fi

# Set MPI wrappers
export MPICC=`which mpicc`
export MPICXX=`which mpicxx`
export MPIF90=`which mpif90`

export USE_HWLOC=OFF
export HWLOC_LIBS=-lhwloc

export HDF5_LIBS="-L${HDF5_ROOT}/lib;-lhdf5_hl;-lhdf5;-lz;-ldl"
#export NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-L${HDF5_ROOT}/lib;-lboost_program_options;-lboost_system;-lnetcdf;-lpnetcdf;-lhdf5_hl;-lhdf5;-lz;-ldl"
#export HDF5_LIBS="-L${HDF5_ROOT}/lib;${HDF5_ROOT}/lib/libhdf5_hl.a;${HDF5_ROOT}/lib/libhdf5.a;-lz;-ldl"
export NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;${BOOST_ROOT}/lib/libboost_program_options.a;${BOOST_ROOT}/lib/libboost_system.a;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_LIBS}"

#export NETCDF_LIBS="-L${BOOST_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${NETCDF_ROOT}/lib;-L${PNETCDF_ROOT}/lib;-lboost_program_options;-lboost_system;${NETCDF_ROOT}/lib/libnetcdf.a;${PNETCDF_ROOT}/lib/libpnetcdf.a;${HDF5_LIBS}"


unset ATTB_ENV

export MPI_POST_FLAG="-map-by;socket:PE=16;--oversubscribe"
