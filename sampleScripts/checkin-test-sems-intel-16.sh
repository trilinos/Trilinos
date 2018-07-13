#!/bin/bash

# This script is meant to be symlinked to the

if [ "$TRILINOS_DIR" == "" ] ; then
  # First try standard directory configurations
  if [ -e ../../Trilinos ] ; then
    TRILINOS_DIR=../../Trilinos
  else
    # Grab from the symlink (only works on Linux)
    _ABS_FILE_PATH=`readlink -f $0` || \
     echo "Could not follow symlink to set TRILINOS_DIR!"
    if [ "$_ABS_FILE_PATH" != "" ] ; then
      _SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
      TRILINOS_DIR=$_SCRIPT_DIR/..
    fi
  fi
fi

echo "TRILINOS_DIR = '$TRILINOS_DIR'"
if [ "$TRILINOS_DIR" == "" ] ; then
  echo "ERROR: Cannot determine TRILINOS_DIR (you must be on a non-Linux system or you must have copied the script instead of symlinking it as per instructions). If you want to try to use this script on this system then please use standard directory structure or set TRILINOS_DIR manually in the env!  But this script and build process is currently only set up to support RHEL Linux 6 machines using the SEMS env."
  exit 1
fi

echo "Loading SEMS modules"
module purge
module load sems-env
module load sems-intel/16.0.3
module load sems-openmpi/1.10.1
module load sems-python/2.7.9
module load sems-cmake/3.5.2
module load sems-git/2.10.1
module load sems-boost/1.63.0/base
module load sems-yaml_cpp/0.5.3/base
module load sems-zlib/1.2.8/base
module load sems-hdf5/1.8.12/parallel
module load sems-netcdf/4.3.2/parallel
module load sems-parmetis/4.0.3/parallel
module load sems-superlu/4.3/base
echo "SEMS modules loaded"

CXX=icpc
CC=icc

######################################################################
# Default configuration options for all builds
######################################################################
rm -f COMMON.config

echo "-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-D BUILD_SHARED_LIBS:BOOL=ON
-D Trilinos_SHOW_DEPRECATED_WARNINGS:BOOL=ON
-D Trilinos_ENABLE_Fortran:BOOL=OFF
-D TPL_ENABLE_yaml-cpp:BOOL=ON
-D TPL_ENABLE_SuperLU:BOOL=ON
-D TPL_ENABLE_Zlib:BOOL=ON
-D TPL_ENABLE_Netcdf:BOOL=ON
-D TPL_ENABLE_HDF5:BOOL=ON
-D TPL_ENABLE_ParMETIS:BOOL=ON
-D TPL_ENABLE_Boost:BOOL=ON
-D TPL_ENABLE_BoostLib:BOOL=ON
-D Trilinos_CXX11_FLAGS:STRING=\"-std=c++11\"
-D CMAKE_CXX_FLAGS:STRING=\"-Wall\"
-D Trilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/sems/SEMSDevEnv.cmake" > COMMON.config

#
# Configuration options for MPI debug build with all complex Scalar
# types disabled
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
" > MPI_RELEASE_DEBUG_INTEL.config

#
# Configuration options for MPI debug build with
# Scalar={double,float,std::complex<double>} enabled, and with only
# the Kokkos::Threads (Kokkos_ENABLE_Pthread) execution space enabled.
#
cat MPI_RELEASE_DEBUG_INTEL.config > MPI_RELEASE_DEBUG_COMPLEX_INTEL.config

echo "
-D Teuchos_ENABLE_FLOAT:BOOL=ON
-D Teuchos_ENABLE_COMPLEX:BOOL=ON
-D Tpetra_INST_FLOAT:BOOL=ON
-D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=ON
" >> MPI_RELEASE_DEBUG_COMPLEX_INTEL.config

_LOCAL_CHECKIN_TEST_DEFAULTS=local-checkin-test-defaults.py
if [ -f $_LOCAL_CHECKIN_TEST_DEFAULTS ] ; then
  echo "File $_LOCAL_CHECKIN_TEST_DEFAULTS already exists, leaving it!"
else
  echo "Creating default file $_LOCAL_CHECKIN_TEST_DEFAULTS!"
  echo "
defaults = [
  \"-j4\",
  \"--ctest-timeout=300\",
  \"--disable-packages=PyTrilinos,Claps,TriKota\",
  \"--skip-case-no-email\",
  \"--default-builds=\",
  \"--extra-builds=MPI_RELEASE_DEBUG_COMPLEX_INTEL\",
  ]
  " > $_LOCAL_CHECKIN_TEST_DEFAULTS
fi

$TRILINOS_DIR/cmake/tribits/ci_support/checkin-test.py "$@"
