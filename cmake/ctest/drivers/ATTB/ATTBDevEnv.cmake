#
# Base options for all GCC 4.8.4 builds on ATTB machines
#

IF ("$ENV{ATTB_ENV}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error: The env var ATTB_ENV must be defined and nonzero"
    " in order to make sure that the ATTB env is set up correctly.")
ENDIF()

MESSAGE("Setting compilers and TPL paths for ATTB system ...")

# Define cmpilers
ASSERT_DEFINED(ENV{MPICC})
ASSERT_DEFINED(ENV{MPICXX})
ASSERT_DEFINED(ENV{MPIF90})
SET(CMAKE_C_COMPILER "$ENV{MPICC}" CACHE FILEPATH
  "Set in gcc-4.8.4-base-options.cmake")
SET(CMAKE_CXX_COMPILER "$ENV{MPICXX}" CACHE FILEPATH
  "Set in gcc-4.8.4-base-options.cmake")
SET(CMAKE_Fortran_COMPILER "$ENV{MPIf90}" CACHE FILEPATH
  "Set in gcc-4.8.4-base-options.cmake")

# Point to the right MPI
ASSERT_DEFINED(ENV{MPI_ROOT})
SET(MPI_BASE_DIR "$ENV{MPI_ROOT}" CACHE PATH
  "Set in gcc-4.8.4-base-options.cmake")

# Turn on explicit template instantaition by default
SET(Trilinos_ENABLE_EXPLICIT_INSTANTIATION  ON  CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Turn of 'int' as a global ordinal
SET(Tpetra_INST_INT_INT OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Turn off float and complex
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Sacado_ENABLE_COMPLEX OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Thyra_ENABLE_COMPLEX OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Tpetra_INST_COMPLEX_DOUBLE OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Tpetra_INST_COMPLEX_FLOAT OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")
SET(Anasazi_ENABLE_COMPLEX OFF CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Enable configure timing
SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING ON CACHE BOOL
  "Set by default in gcc-4.8.3-base-options.cmake")

# Set up TPL stuff

# BLAS
#ASSERT_DEFINED(ENV{BLAS_ROOT})
#SET(BLAS_LIBRARY_DIRS "$ENV{BLAS_ROOT}/lib"
#  CACHE PATH "Set in ATTBDevEnv.cmake")

# LAPACK
#ASSERT_DEFINED(ENV{LAPACK_ROOT})
#SET(LAPACK_LIBRARY_DIRS "$ENV{LAPACK_ROOT}/lib"
#  CACHE PATH "Set in ATTBDevEnv.cmake")

# Boost
#ASSERT_DEFINED(ENV{BOOST_ROOT})
#SET(Boost_INCLUDE_DIRS "$ENV{Boost_ROOT}/include"
#  CACHE PATH "Set in ATTBDevEnv.cmake")
#SET(Boost_LIBRARY_DIRS "$ENV{Boost_ROOT}/lib"
#  CACHE PATH "Set in ATTBDevEnv.cmake")

# BoostLib
#ASSERT_DEFINED(ENV{BOOST_ROOT})
#SET(BoostLib_INCLUDE_DIRS "$ENV{Boost_ROOT}/include"
#  CACHE PATH "Set in ATTBDevEnv.cmake")
#SET(BoostLib_LIBRARY_DIRS "$ENV{Boost_ROOT}/lib"
#  CACHE PATH "Set in ATTBDevEnv.cmake")

# Netcdf
#ASSERT_DEFINED(ENV{NETCDF_ROOT})
#ASSERT_DEFINED(ENV{HDF5_ROOT})
#SET(TPL_Netcdf_INCLUDE_DIRS "$ENV{NETCDF_ROOT}/include;$ENV{HDF5_ROOT}/include"
#  CACHE PATH "Set in ATTBDevEnv.cmake"  CACHE PATH "Set in ATTBDevEnv.cmake")
#SET(Netcdf_LIBRARY_DIRS "$ENV{NETCDF_ROOT}/lib;$ENV{HDF5_ROOT}/lib"
#  CACHE PATH "Set in ATTBDevEnv.cmake")
#SET(Netcdf_LIBRARY_NAMES "netcdf;hdf5_hl;hdf5_fortran;hdf5;z"
#  CACHE PATH "Set in ATTBDevEnv.cmake")
