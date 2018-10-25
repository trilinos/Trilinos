#!/bin/bash
#
# Trilinos check-in test script that uses the SEMS modules
# Time-stamp: <2016-09-22 15:56:59 mhoemme>
#
# DO NOT use this file without first reading it and editing some parts
# as necessary!
#
# I. This script uses SEMS modules
#
# This bash script wraps Trilinos' check-in test script.  It uses
# Sandia's SEMS "module" system to load compilers and third-party
# libraries that Trilinos' build system uses.  This script will NOT
# work with systems that do not use the same module system.  If your
# computer has some kind of module-based system, it would likely not
# be hard to modify this script to make it work.  If you do not have a
# module-based system at all, this script is likely not for you.
#
# This script writes to a file called "modules-used-by-checkin-tests"
# in this directory.  That file gets a list of all the module-related
# commands that this script runs.  You may run the contents of
# modules-used-by-checkin-tests as a Bash script, in order to load up
# the same modules.  Please note that both this script and
# modules-used-by-checkin-tests first run "modules purge", which
# clears the currently loaded set of modules.
#
# II. How to run this script
# II.A. Run this script by "source"ing it
#
# In order to pick up the SEMS modules correctly, you should always
# run this script by "sourcing" it, like this (replacing the ellipsis
# with any command-line arguments):
#
#   $ source ./checkin-test-using-sems-modules.sh ...
#
# II.B. Command-line arguments
# II.B.1. If you changed something in Trilinos
#
# If you change something in some package, have made commits, and
# would like to push to Trilinos, do the following (replace $N with a
# positive integer indicating desired maximum build parallelism):
#
#   $ source ./checkin-test-using-sems-modules.sh --do-all --push \
#       --make-options="-j$N"
#
# If you are not ready to push and not confident about commits, try
# this instead:
#
#   $ source ./checkin-test-using-sems-modules.sh --configure \
#       --build --test --allow-no-pull --enable-all-packages=off \
#       --make-options="-j$N"
#
# The options "--configure --build --test" enable all the things that
# Trilinos needs to do for a full test.  The "--allow-no-pull" option
# convinces Trilinos not to attempt a pull or merge.  (This is a must
# if you don't have a network connection.)  The
# "--enable-all-packages=off" option does not actually disable
# packages or prevent them from being enabled.  It just convinces the
# check-in test script not to worry about the "--allow-no-pull"
# option.
#
# You may also wish to test packages that Trilinos does not enable by
# default.  In order to do that, add the following to the above
# command:
#
#   --default-builds="" --extra-builds=MPI_DEBUG,SERIAL_RELEASE
#
# All of the above commands just test the MPI_DEBUG and SERIAL_RELEASE
# builds (see below; "SERIAL" in this context just means "no MPI," and
# has nothing to do with thread parallelism).  If you want to test
# other builds, add them to the --extra-builds argument.  For example:
#
#   --default-builds="" --extra-builds=MPI_RELEASE,SERIAL_DEBUG
#
# Those builds _must_ exist in the check-in script (they do in this
# case; see below).  You may create as many builds as you have disk
# space to store.  See below in this script for examples.  We set
# --default-builds="" because, by default, --default-builds always
# includes MPI_DEBUG and SERIAL_RELEASE.
#
# II.B.2. If you have not yet changed anything in Trilinos
#
# Trilinos enables packages based on what files you changed.  If you
# have not changed any files, Trilinos will not enable any packages by
# default.  In that case, if you want to enable one or more packages,
# you must enable them explicitly, as a comma-delimited list given to
# the "--enable-packages" command-line option.  For example:
# "--enable-packages=Tpetra,Ifpack2,Amesos2".  Package names are case
# sensitive.  Trilinos/PackagesList.cmake lists all packages.
#
# You may also include subpackages in the list you give to the
# --enable-packages option.  For a list of all valid subpackage names,
# look at the output of Trilinos' configuration process (e.g.,
# MPI_DEBUG/configure.out).  What follows after "Final set of
# non-enabled SE packages" shows all non-enabled subpackages and
# packages, and what follows after "Final set of enabled SE packages"
# shows all enabled subpackages and packages.
#
# II.C. How to get help
#
# For detailed documentation on all command-line arguments, just run
# this script with the "--help" command-line option, like this:
#
#   $ source .checkin-test-using-sems-modules.sh --help
#

# Pick up command-line arguments (DO NOT CHANGE THIS)
EXTRA_ARGS=$@

# This is where your Trilinos source check-out lives.
TRILINOS_PATH=/scratch/prj/Trilinos/Trilinos

# Where the SEMS modules live.
#
# We use a "modules" system to load compilers and libraries.  SEMS
# (some acronym refering to our software infrastructure team)
# maintains a modules system.  Your computer must mount this directory
# in order to load modules from it.
TPL_PATH=/projects/install/rhel6-x86_64/sems/tpl

# Whether to use ninja, instead of make, for builds.  The ninja binary
# must be in your PATH if you enable this option.
#
# NOTE: This currently ONLY works if you run ninja by hand.  That is,
# it only works if you only use the check-in test script with
# --configure.  If you use --build (or --do-all or --local-do-all), it
# won't work.  You'll have to use make in that case.

#USE_NINJA=ON
USE_NINJA=OFF

#
# Decisions about compilers
#

# Host compiler chain (clang, gcc, or intel).  This must match the
# host compiler version (see below).

HOST_COMPILER_CHAIN=clang
#HOST_COMPILER_CHAIN=gcc
#HOST_COMPILER_CHAIN=intel

# Host compiler version.  This must match the host compiler chain (see
# above).

HOST_COMPILER_VERSION=3.6.1
#HOST_COMPILER_VERSION=4.7.2
#HOST_COMPILER_VERSION=5.3.0
#HOST_COMPILER_VERSION=16.0.3

# If you want to enable CUDA, set CUDA=ON right now.

#CUDA=ON
CUDA=OFF

# GCC and Intel both implement OpenMP.  Only Clang >= 3.8 implements
# OpenMP.  You may enable both OpenMP and CUDA at the same time.  I
# generally always enable OpenMP as long as the compiler supports it.
# This is because many Trilinos users already enable OpenMP in their
# builds.  Compilers may behave differently when OpenMP is enabled,
# even if the code doesn't use OpenMP.
if [ "${HOST_COMPILER_CHAIN}" == "clang" ]; then
  OPENMP=OFF
else
  OPENMP=ON
fi

# MPI implementation and version
#
# If you enable CUDA, this script requires "openmpi" (OpenMPI) as the
# MPI implementation.  This is because, when CUDA is enabled, this
# script sets the OMPI_CXX environment variable in order to use
# "nvcc_wrapper" (see explanation below) as the "compiler" that MPI's
# compiler wrappers invoke.

MPI_IMPLEMENTATION=openmpi
#MPI_VERSION=1.8.7
MPI_VERSION=1.10.1

# The CUDA version always needs to be defined, even if not using CUDA.
# This is because the SEMS module system always wants to load CUDA
# when loading a compiler.
CUDA_VERSION=7.5.18

# Trilinos requires at least CMake 2.8.11.  3.x works fine.  In fact,
# we recommend 3.3.2 if it is available, for best performance.
CMAKE_VERSION=3.3.2

# You don't need to set any of these.
CMAKE_MODULE=sems-cmake/${CMAKE_VERSION}
COMPILER_SUFFIX=${HOST_COMPILER_CHAIN}/${HOST_COMPILER_VERSION}
COMPILER_MODULE=sems-${COMPILER_SUFFIX}
CUDA_SUFFIX=cuda/${CUDA_VERSION}
CUDA_MODULE=${CUDA_SUFFIX}
MPI_SUFFIX=${MPI_IMPLEMENTATION}/${MPI_VERSION}
MPI_MODULE=sems-${MPI_SUFFIX}

#
# Decisions about third-party libraries (TPLs)
#
# Trilinos depends on TPLs.  Almost all dependencies are optional.
# However, some optional Trilinos features may depend on certain TPLs.
# Generally, the more TPLs you enable, the more Trilinos features
# you'll enable, and thus, the more Trilinos code you'll build and
# test.
#
# I've added variables HAVE_$TPL for all values of $TPL.  You may set
# any of these to OFF if you want to disable that TPL.  It's OK to set
# $TPL_LIBRARY_DIRS, etc., even if TPL_ENABLE_$TPL is OFF; the latter
# overrides the former.
#

# Boost implements many different C++ utilities.  Trilinos separates
# Boost into two TPLs: the header-only part of Boost ("Boost" TPL),
# and the part of Boost that also requires linking against a library
# ("BoostLib" TPL).
HAVE_BOOST=ON
HAVE_BOOSTLIB=ON
BOOST_VERSION=1.59.0
BOOST_SUFFIX=boost/${BOOST_VERSION}
BOOST_MODULE=sems-${BOOST_SUFFIX}
BOOST_PATH=${TPL_PATH}/${BOOST_SUFFIX}/${COMPILER_SUFFIX}/base

# SuperLU provides a single-MPI-process sparse LU factorization.
HAVE_SUPERLU=ON
SUPERLU_VERSION=4.3
SUPERLU_SUFFIX=superlu/${SUPERLU_VERSION}
SUPERLU_MODULE=sems-${SUPERLU_SUFFIX}
SUPERLU_PATH=${TPL_PATH}/${SUPERLU_SUFFIX}/${COMPILER_SUFFIX}/base

# NetCDF implements parallel file I/O.
HAVE_NETCDF=ON
NETCDF_VERSION=4.3.2
NETCDF_SUFFIX=netcdf/${NETCDF_VERSION}
NETCDF_MODULE=sems-${NETCDF_SUFFIX}
NETCDF_PATH=${TPL_PATH}/${NETCDF_SUFFIX}/${COMPILER_SUFFIX}/${MPI_SUFFIX}

# HDF5 implements parallel file I/O.
HAVE_HDF5=ON
HDF5_VERSION=1.8.12
HDF5_SUFFIX=hdf5/${HDF5_VERSION}
HDF5_MODULE=sems-${HDF5_SUFFIX}
HDF5_PATH=${TPL_PATH}/${HDF5_SUFFIX}/${COMPILER_SUFFIX}/${MPI_SUFFIX}

# zlib implements data compression.
HAVE_ZLIB=OFF
ZLIB_VERSION=1.2.8
ZLIB_SUFFIX=zlib/${ZLIB_VERSION}
ZLIB_MODULE=sems-${ZLIB_SUFFIX}
ZLIB_PATH=${TPL_PATH}/${ZLIB_SUFFIX}/${COMPILER_SUFFIX}/base

# ParMETIS implements parallel graph partitioning.
HAVE_PARMETIS=ON
PARMETIS_VERSION=4.0.3
PARMETIS_SUFFIX=parmetis/${PARMETIS_VERSION}
PARMETIS_MODULE=sems-${PARMETIS_SUFFIX}
PARMETIS_PATH=${TPL_PATH}/${PARMETIS_SUFFIX}/${COMPILER_SUFFIX}/${MPI_SUFFIX}

#
# Load the modules
#

# Last error checking before we change the environment.
if [ "${CUDA}" == 'ON' -a "${MPI_IMPLEMENTATION}" != 'openmpi' ]; then
  echo "If you enable CUDA, this script only works with OpenMPI."
fi
if [ "${HOST_COMPILER_CHAIN}" != "gcc" -a "${HOST_COMPILER_CHAIN}" != "intel" -a "${HOST_COMPILER_CHAIN}" != "clang" ]; then
  echo "I don't know about the host compiler called \"${HOST_COMPILER_CHAIN}\"."
fi

# We start by clearing any loaded modules.
module purge

# The SEMS module system requires that this module be loaded.
module load sems-env

module load ${CMAKE_MODULE}
module load ${COMPILER_MODULE}
module load ${MPI_MODULE}

if [ "${CUDA}" == "ON" ]; then
  # Load modules needed for CUDA builds.
  module load ${CUDA_MODULE}
  # Set environment variables that Trilinos needs for building and
  # running with CUDA.  "nvcc_wrapper" is a script that makes NVCC
  # behave more like the other compilers.  Note that OMPI_CXX only
  # works if you are using OpenMPI.
  export OMPI_CXX=${TRILINOS_PATH}/packages/kokkos/config/nvcc_wrapper
  export CUDA_LAUNCH_BLOCKING=1
  export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1
fi

# Load optional third-party libraries (TPLs).
module load ${BOOST_MODULE}
module load ${SUPERLU_MODULE}
module load ${NETCDF_MODULE}
module load ${HDF5_MODULE}
module load ${ZLIB_MODULE}
module load ${PARMETIS_MODULE}

#
# Generate the "modules-used-by-checkin-tests" file
#

# The first line doesn't append, so that we clear out the file.
echo "module purge" > modules-used-by-checkin-tests
# Subsequent lines append.
echo "module load sems-env" >> modules-used-by-checkin-tests
echo "module load ${CMAKE_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${COMPILER_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${MPI_MODULE}" >> modules-used-by-checkin-tests
if [ "${CUDA}" == "ON" ]; then
  echo "module load ${CUDA_MODULE}" >> modules-used-by-checkin-tests
  echo "export OMPI_CXX=${TRILINOS_PATH}/packages/kokkos/config/nvcc_wrapper" >> modules-used-by-checkin-tests
  echo "export CUDA_LAUNCH_BLOCKING=1" >> modules-used-by-checkin-tests
  echo "export CUDA_MANAGED_FORCE_DEVICE_ALLOC=1" >> modules-used-by-checkin-tests
fi

echo "module load ${BOOST_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${SUPERLU_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${NETCDF_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${HDF5_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${ZLIB_MODULE}" >> modules-used-by-checkin-tests
echo "module load ${PARMETIS_MODULE}" >> modules-used-by-checkin-tests

# "Serial" (non-MPI) builds need the CXX and CC variables.  If you
# need to support a compiler other than GCC, Intel, or Clang, you'll
# need to set these by hand.
if [ "${HOST_COMPILER_CHAIN}" == "gcc" ]; then
  CXX=g++
  CC=gcc
elif [ "${HOST_COMPILER_CHAIN}" == "intel" ]; then
  CXX=icpc
  CC=icc
elif [ "${HOST_COMPILER_CHAIN}" == "clang" ]; then
  CXX=clang++
  CC=clang
else
  echo "I don't know about the host compiler called \"${HOST_COMPILER_CHAIN}\"."
  exit -1
fi

######################################################################
# Set up options for builds
######################################################################

# Trilinos' check-in test script can test multiple configurations
# ("builds") of Trilinos.  By default, it tests two configurations:
# MPI_DEBUG (MPI enabled, debug build) and SERIAL_RELEASE (MPI
# disabled, release build).  You may add other builds to that list.
#
# You, the user of the check-in test script, are responsible for
# defining how Trilinos should set up those builds.  This includes
# what compilers to use, what third-party libraries (TPL) to link in
# and where to find them, and what Trilinos options to set.  Trilinos
# takes those options for each build from different files in this
# directory.  It first looks in a file called COMMON.config for the
# default options for _all_ builds.  Then, for each build $BUILD, it
# looks in a file called $BUILD.config for that build's options.  This
# script creates those files (by writing to them using "echo"), so
# that Trilinos' check-in test script can find and use them.

######################################################################
# Default configuration options for all builds
######################################################################

if [ "${USE_NINJA}" == 'ON' ]; then
  echo "-G Ninja" > COMMON.config
else
  # Do this so we can append to COMMON.config below.
  rm -f COMMON.config
fi

# Trilinos_ENABLE_EXPLICIT_INSTANTIATION: "Explict instantiation"
# refers to C++ explicit template instantiation (ETI).  Enabling this
# option almost always improves build time.
#
# BUILD_SHARED_LIBS: ON tells Trilinos to build dynamic shared
# libraries (.so, .dylib, etc.).  OFF tells Trilinos to build static
# libraries.  Dynamic shared libraries improve build time and use less
# disk space.
#
# CMAKE_CXX_FLAGS: (Extra) flags for the C++ compiler.  Trilinos
# considers warnings bugs, so we set -Wall --pedantic.  You may need
# to change these flags if you change compilers.  We demonstrate this
# in the if/else/fi statement below that appends to COMMON.config.
#
# Trilinos_ENABLE_OpenMP: Setting option takes care of adding any
# OpenMP-related flags to the compiler's or linker's command line.
# You generally do not need to add such flags to CMAKE_CXX_FLAGS.
#
# Trilinos_SHOW_DEPRECATED_WARNINGS: Trilinos lets developers mark
# functions, classes, etc. as "deprecated."  If the compiler supports
# this option, then it will emit deprecated warnings whenever some
# code uses the deprecated thing.
#
# Trilinos_ENABLE_Fortran: Trilinos as a whole does not require
# Fortran.  Some packages may require Fortran, and others may have
# optional features that require Fortran.  Clang does not currently
# provide a Fortran compiler, so for the sake of uniformity across
# compiler options in this script, we disable Fortran here.
#
# Teuchos_ENABLE_COMPLEX and Tpetra_INST_COMPLEX_DOUBLE: If you want
# to enable Scalar=std::complex<double> in Tpetra, you must set these
# two options.  For more Tpetra configure-time options, see
# Trilinos/packages/tpetra/doc/FAQ.txt.
#
# TPL-related options: TPL_ENABLE_${TPL} enables (or disables) the
# TPL.  Trilinos must know about the TPL in order for this to work.
# (For the list of all supported TPLs, see Trilinos/TPLsList.cmake.)
# If you just enable a TPL without telling Trilinos where to look for
# its include files or libraries, then Trilinos will try to guess
# where the TPL's things might live.  If you use the SEMS modules, you
# need to tell Trilinos explicitly where those include files and
# libraries live, by setting the ${TPL}_INCLUDE_DIRS and
# ${TPL}_LIBRARY_DIRS options.
#
# It may make sense to have specific builds that enable or disable
# certain TPLs, or that use different versions of the same TPLs.  This
# may motivate moving some of the TPL-related options out of
# COMMON.config, and into specific builds' $BUILD.config files.
#
echo "
-D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
-D BUILD_SHARED_LIBS:BOOL=ON
-D TPL_ENABLE_CUDA=${CUDA}
  -D Kokkos_ENABLE_Cuda_UVM:BOOL=${CUDA}
  -D Tpetra_INST_CUDA=${CUDA}
-D Trilinos_ENABLE_OpenMP:BOOL=${OPENMP}
-D Trilinos_SHOW_DEPRECATED_WARNINGS:BOOL=ON
-D Trilinos_ENABLE_Fortran:BOOL=OFF
-D Teuchos_ENABLE_COMPLEX:BOOL=ON
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=ON
-D TPL_ENABLE_SuperLU:BOOL=${HAVE_SUPERLU}
  -D SuperLU_INCLUDE_DIRS:FILEPATH=\"${SUPERLU_PATH}/include\"
  -D SuperLU_LIBRARY_DIRS:FILEPATH=\"${SUPERLU_PATH}/lib\"
-D TPL_ENABLE_Zlib:BOOL=${HAVE_ZLIB}
  -D Zlib_INCLUDE_DIRS:FILEPATH=\"${ZLIB_PATH}/include\"
  -D Zlib_LIBRARY_DIRS:FILEPATH=\"${ZLIB_PATH}/lib\"
-D TPL_ENABLE_Netcdf:BOOL=${HAVE_NETCDF}
  -D Netcdf_INCLUDE_DIRS:FILEPATH=\"${NETCDF_PATH}/include\"
  -D Netcdf_LIBRARY_DIRS:FILEPATH=\"${NETCDF_PATH}/lib\"
-D TPL_ENABLE_HDF5:BOOL=${HAVE_HDF5}
  -D HDF5_INCLUDE_DIRS:FILEPATH=\"${HDF5_PATH}/include\"
  -D HDF5_LIBRARY_DIRS:FILEPATH=\"${HDF5_PATH}/lib\"
-D TPL_ENABLE_ParMETIS:BOOL=${HAVE_PARMETIS}
  -D ParMETIS_INCLUDE_DIRS:FILEPATH=\"${PARMETIS_PATH}/include\"
  -D ParMETIS_LIBRARY_DIRS:FILEPATH=\"${PARMETIS_PATH}/lib\"
-D TPL_ENABLE_Boost:BOOL=${HAVE_BOOST}
  -D Boost_INCLUDE_DIRS=${BOOST_PATH}/include
  -D Boost_LIBRARY_DIRS=${BOOST_PATH}/lib
-D TPL_ENABLE_BoostLib:BOOL=${HAVE_BOOSTLIB}
  -D BoostLib_INCLUDE_DIRS=${BOOST_PATH}/lib
  -D BoostLib_LIBRARY_DIRS=${BOOST_PATH}/lib
" >> COMMON.config

# ^^^ We append (>>) rather than overwrite (>) for COMMON.config only,
# because we want the CMake generator option (-G Ninja) to come before
# other arguments passed to CMake.  (At least it makes the resulting
# CMake command easier to read.)

#
# Compiler-specific options to add to the list of default
# configuration options for all builds.
#

if [ "${CUDA}" == "ON" ]; then
  echo "-D CMAKE_CXX_FLAGS:STRING=\"-Wall\"" >> COMMON.config
elif [ "${HOST_COMPILER_CHAIN}" == "intel" ]; then
  # Intel (<= 15?) needs a special flag to enable C++11.
  echo "-D Trilinos_CXX11_FLAGS:STRING=\"-std=c++11\"" >> COMMON.config
  # If you try to give the Intel compiler the --pedantic command-line
  # option, it may report a warning like the following:
  #
  # icpc: command line warning #10006: ignoring unknown option '-fpedantic'
  #
  # As a result, we only include this flag when _not_ using the Intel compiler.
  echo "-D CMAKE_CXX_FLAGS:STRING=\"-Wall\"" >> COMMON.config
else
  # Both GCC and Clang, but not Intel, understand the "--pedantic" option.
  echo "-D CMAKE_CXX_FLAGS:STRING=\"-Wall --pedantic\"" >> COMMON.config
fi

#
# Configuration options for MPI debug build
#
# CMAKE_BUILD_TYPE: "DEBUG" or "RELEASE".  In this case, we want
# "DEBUG".
#
# TPL_ENABLE_MPI: Trilinos treats MPI sort of like a TPL.  If Trilinos
# can find your MPI implementation easily, then you don't have to
# specify more than this.  Otherwise, you may need to give it more
# information.  See other scripts for examples.
#
# Kokkos_ENABLE_DEBUG: If set, turn on bounds checking for
# Kokkos::View objects, as well as other expensive but useful run-time
# error checks relating to Kokkos.
#
# Teuchos_ENABLE_DEBUG: If set, turn on bounds checking for
# Teuchos::Array* objects, as well as other expensive but useful
# run-time error checks relating to Teuchos.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
" > MPI_DEBUG.config

#
# Configuration options for MPI debug build with ETI OFF.
#
# Turning off ETI usually increases build time, so disable complex and
# bounds-checking ("debug") options to keep build time down.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=OFF
-D Teuchos_ENABLE_DEBUG:BOOL=OFF
-D Teuchos_ENABLE_COMPLEX:BOOL=OFF
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=OFF
" > MPI_DEBUG_NO_ETI.config

#
# Configuration options for MPI debug build with all complex Scalar
# types disabled
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_COMPLEX:BOOL=OFF
-D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=OFF
" > MPI_DEBUG_REAL.config

#
# Configuration options for MPI debug build with GO=int OFF
#
# This may need the following:
#   -D Amesos2_ENABLE_Epetra:BOOL=OFF
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
-D Tpetra_INST_INT_INT:BOOL=OFF
" > MPI_DEBUG_NO_INT.config

#
# Configuration options for MPI debug build with static libraries
#
# Setting BUILD_SHARED_LIBS=OFF here overrides the default, which we
# set in COMMON.config (see above).
#
echo "
-D BUILD_SHARED_LIBS:BOOL=OFF
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
" > MPI_DEBUG_STATIC.config

#
# Configuration options for MPI debug build with
# Scalar={double,float,std::complex<double>} enabled, and with only
# the Kokkos::Threads (Kokkos_ENABLE_Pthread) execution space enabled.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
-D Kokkos_ENABLE_Serial:BOOL=OFF
-D Trilinos_ENABLE_OpenMP:BOOL=OFF
-D Kokkos_ENABLE_OpenMP:BOOL=OFF
-D Kokkos_ENABLE_Pthread:BOOL=ON
-D Teuchos_ENABLE_COMPLEX:BOOL=ON
-D Tpetra_INST_FLOAT:BOOL=ON
-D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=ON
" > MPI_DEBUG_COMPLEX.config

#
# Configuration options for MPI debug build (secondary stable included)
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
" > MPI_DEBUG_SS.config

#
# Configuration options for MPI release build.
#
# This build disables use of OpenMP in Tpetra, and disables GO=int,
# all to save build time.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D TPL_ENABLE_MPI:BOOL=ON
-D Tpetra_INST_OPENMP:BOOL=OFF
-D Tpetra_INST_SERIAL:BOOL=ON
-D Tpetra_INST_INT_INT:BOOL=OFF
" > MPI_RELEASE.config

#
# Configuration options for MPI release build where Kokkos::Serial is
# the default execution space, and GO=long is enabled in Tpetra.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D TPL_ENABLE_MPI:BOOL=ON
-D Kokkos_ENABLE_Serial:BOOL=ON
-D Tpetra_INST_SERIAL:BOOL=ON
-D Trilinos_ENABLE_OpenMP:BOOL=OFF
-D Kokkos_ENABLE_OpenMP:BOOL=OFF
-D Tpetra_INST_OPENMP:BOOL=OFF
-D Kokkos_ENABLE_Pthread:BOOL=OFF
-D Tpetra_INST_PTHREAD:BOOL=OFF
-D Tpetra_INST_FLOAT:BOOL=OFF
-D Tpetra_INST_COMPLEX_FLOAT:BOOL=OFF
-D Tpetra_INST_COMPLEX_DOUBLE:BOOL=OFF
-D Tpetra_INST_INT_LONG:BOOL=ON
" > MPI_RELEASE_SERIAL_LONG.config

#
# Configuration options for MPI release build (secondary stable included)
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D TPL_ENABLE_MPI:BOOL=ON
" > MPI_RELEASE_SS.config

#
# Configuration options for serial release build
#
# Disable HDF5 and ParMETIS, since the SEMS installation of that TPL
# seems to require MPI.  If ParMETIS is enabled, code that uses Zoltan
# fails to link.  If HDF5 is enabled, EpetraExt's HDF5 interface fails
# to build.  Perhaps I didn't pick up the right modules.
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D CMAKE_CXX_COMPILER:FILEPATH=\"${CXX}\"
-D CMAKE_C_COMPILER:FILEPATH=\"${CC}\"
-D TPL_ENABLE_HDF5=OFF
-D TPL_ENABLE_ParMETIS=OFF
" > SERIAL_RELEASE.config

#-D TPL_ENABLE_Pthread:BOOL=ON
#-D Kokkos_ENABLE_Pthread:BOOL=OFF

#
# Configuration options for serial release build with static libraries
#
echo "
-D BUILD_SHARED_LIBS:BOOL=OFF
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D CMAKE_CXX_COMPILER:FILEPATH=\"${CXX}\"
-D CMAKE_C_COMPILER:FILEPATH=\"${CC}\"
-D TPL_ENABLE_HDF5=OFF
-D TPL_ENABLE_ParMETIS=OFF
" > SERIAL_RELEASE_STATIC.config

#
# Configuration options for serial release build (secondary stable included)
#
echo "
-D CMAKE_BUILD_TYPE:STRING=RELEASE
-D CMAKE_CXX_COMPILER:FILEPATH=\"${CXX}\"
-D CMAKE_C_COMPILER:FILEPATH=\"${CC}\"
-D TPL_ENABLE_HDF5=OFF
-D TPL_ENABLE_ParMETIS=OFF
" > SERIAL_RELEASE_SS.config

#
# Configuration options for serial debug build
#
echo "
-D CMAKE_BUILD_TYPE:STRING=DEBUG
-D CMAKE_CXX_COMPILER:FILEPATH=\"${CXX}\"
-D CMAKE_C_COMPILER:FILEPATH=\"${CC}\"
-D Kokkos_ENABLE_DEBUG:BOOL=ON
-D Teuchos_ENABLE_DEBUG:BOOL=ON
-D Kokkos_ENABLE_Serial:BOOL=OFF
-D Kokkos_ENABLE_OpenMP:BOOL=OFF
-D Kokkos_ENABLE_Pthread:BOOL=ON
-D TPL_ENABLE_HDF5=OFF
-D TPL_ENABLE_ParMETIS=OFF
" > SERIAL_DEBUG.config

# This is where we actually invoke Trilinos' check-in test script.
#
# Any command-line arguments given to the outer (our) script get
# passed along to this (inner) script.  If you have arguments that you
# always want to pass in, you may put them here.  For example,
# --ctest-timeout=$N sets the test timeout (the amount of time a test
# may run before it is stopped and considered failed) to $N seconds.
# The default is 180, which is not enough for some tests (esp. MueLu's
# unit tests).
#
# Here are some other options you may find useful:
#
# --disable-packages=FEI,PyTrilinos,Moertel,STK,SEACAS,OptiPack,Rythmos,Intrepid,ROL,Panzer
#
#   This disables the given list of packages.  You may find this
#   helpful if some packages are known to give you trouble.
#
# --ctest-options="-E '(TeuchosCore_MemoryManagement_UnitTests|KokkosCore_|KokkosContainers_|Tpetra_MultiPrecExample_double_float|Tpetra_RTIInlineCG|TeuchosParameterList_ObjectBuilder_UnitTests)'"
#
#   This disables specific tests, based on pattern-matching.
#   --ctest-options="-E SomeTest" disables all tests that have
#   "SomeTest" (sans quotes) in their name.  Vertical bars (|)
#   separate OR (disjunctive) matches, and '()' encloses a list of
#   such disjunctions.

${TRILINOS_PATH}/checkin-test.py \
--ctest-timeout=400 \
$EXTRA_ARGS
