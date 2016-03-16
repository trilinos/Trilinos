#
# Base options for all GCC 4.8.4 builds on ATTB machines
#

IF ("$ENV{ATTB_ENV}" STREQUAL "")
  MESSAGE(FATAL_ERROR "Error: The env var ATTB_ENV must be defined and nonzero"
    " in order to make sure that the ATTB env is set up correctly.")
ENDIF()

MESSAGE("-- ATTB: Setting compilers and TPL paths for ATTB system ...")

ASSERT_DEFINED(ENV{GCC_VERSION})
SET(GCC_LIBRARY_PATH
  "/home/projects/x86-64-haswell/gnu/$ENV{GCC_VERSION}/lib/gcc/x86_64-unknown-linux-gnu/$ENV{GCC_VERSION}")
MESSAGE("-- ATTB: Set env var LIBRARY_PATH=.../x86_64-unknown-linux-gnu/$ENV{GCC_VERSION} to avoid problems with"
  " CMake not adding lib paths for other TPLs ...")
SET(ENV{LIBRARY_PATH} "${GCC_LIBRARY_PATH}")
MESSAGE("-- ENV{LIBRARY_PATH} = $ENV{LIBRARY_PATH}")

# Define cmpilers
ASSERT_DEFINED(ENV{MPICC})
ASSERT_DEFINED(ENV{MPICXX})
ASSERT_DEFINED(ENV{MPIF90})
SET(CMAKE_C_COMPILER "$ENV{MPICC}" CACHE FILEPATH
 "Set in ATTBDevEnv.cmake")
SET(CMAKE_CXX_COMPILER "$ENV{MPICXX}" CACHE FILEPATH
 "Set in ATTBDevEnv.cmake")
SET(CMAKE_Fortran_COMPILER "$ENV{MPIf90}" CACHE FILEPATH
 "Set in ATTBDevEnv.cmake")

# Add rpath for compiler libraries and gomp for parts built with OpenMP
IF (TPL_FIND_SHARED_LIBS)
  SET(LDL_LINK_ARG " -lldl")
ELSE()
  SET(LDL_LINK_ARG)
ENDIF()
ASSERT_DEFINED(ENV{GCC_PATH})
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS
  "-lgomp -lgfortran${LDL_LINK_ARG} -ldl"
  CACHE STRING "Set in ATTBDevEnv.cmake")

#  -Wl,-rpath,$ENV{GCC_PATH}/lib64

# Point to the right MPI
ASSERT_DEFINED(ENV{MPI_ROOT})
SET(MPI_BASE_DIR "$ENV{MPI_ROOT}" CACHE PATH
 "Set in ATTBDevEnv.cmake")

# Turn on explicit template instantaition by default
SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION  ON  CACHE BOOL
  "Set in ATTBDevEnv.cmake")

# Turn of 'int' as a global ordinal
#SET(Tpetra_INST_INT_INT OFF CACHE BOOL
#  "Set in ATTBDevEnv.cmake")
# NOT: We have to keep int as a global ordinal in order to support Epetra
# adapters.

# Turn off float and complex
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Sacado_ENABLE_COMPLEX OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Thyra_ENABLE_COMPLEX OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_DOUBLE OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_FLOAT OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(Anasazi_ENABLE_COMPLEX OFF CACHE BOOL "Set in ATTBDevEnv.cmake")

# Enable configure timing
SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING ON CACHE BOOL
  "Set in ATTBDevEnv.cmake")

#
# Set up TPL stuff
#
# These ATTB systems are currently set up to only have static libs.
#
# We disable TPLs that we know are not on this system so that we disable
# downstream SE packages that have required dependencies on these TPLs
# (e.g. some SEACAS and STK subpackages).
#
# We default enable all of the TPLs that we know are on this system.  In this
# way, support for these TPLs will be turned on by default.
#

# Always find and use static libs on this system
SET(TPL_FIND_SHARED_LIBS OFF CACHE BOOL
  "Set in ATTBDevEnv.cmake")
SET(${PROJECT_NAME}_LINK_SEARCH_START_STATIC ON  CACHE BOOL
  "Set in ATTBDevEnv.cmake")

# Disable a bunch of TPLs that are not on this system
SET(TPL_ENABLE_GLM OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(TPL_ENABLE_Matio OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(TPL_ENABLE_SuperLU OFF CACHE BOOL "Set in ATTBDevEnv.cmake")
SET(TPL_ENABLE_X11 OFF CACHE BOOL "Set in ATTBDevEnv.cmake")

# BLAS
SET(TPL_ENABLE_BLAS ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{BLAS_ROOT})
SET(BLAS_LIBRARY_DIRS "$ENV{BLAS_ROOT}/lib"
  CACHE PATH "Set in ATTBDevEnv.cmake")

# LAPACK
SET(TPL_ENABLE_LAPACK ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{LAPACK_ROOT})
SET(LAPACK_LIBRARY_DIRS "$ENV{LAPACK_ROOT}/lib"
  CACHE PATH "Set in ATTBDevEnv.cmake")

# Boost
SET(TPL_ENABLE_Boost ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{BOOST_ROOT})
SET(Boost_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(Boost_LIBRARY_DIRS "$ENV{BOOST_ROOT}/lib"
  CACHE PATH "Set in ATTBDevEnv.cmake")

# BoostLib
SET(TPL_ENABLE_BoostLib ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{BOOST_ROOT})
SET(BoostLib_INCLUDE_DIRS "$ENV{BOOST_ROOT}/include"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(BoostLib_LIBRARY_DIRS "$ENV{BOOST_ROOT}/lib"
  CACHE PATH "Set in ATTBDevEnv.cmake")

# HDF5
SET(TPL_ENABLE_HDF5 ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{HDF5_ROOT})
ASSERT_DEFINED(ENV{ZLIB_ROOT})
SET(TPL_HDF5_INCLUDE_DIRS "$ENV{HDF5_ROOT}/include;$ENV{ZLIB_ROOT}/include"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(HDF5_LIBRARY_DIRS "$ENV{HDF5_ROOT}/lib;$ENV{ZLIB_ROOT}/lib"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5_fortran;hdf5;z"
  CACHE STRING "Set in ATTBDevEnv.cmake")

# Netcdf
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "Set in ATTBDevEnv.cmake")
ASSERT_DEFINED(ENV{NETCDF_ROOT})
ASSERT_DEFINED(ENV{PNETCDF_ROOT})
ASSERT_DEFINED(ENV{HDF5_ROOT})
SET(TPL_Netcdf_INCLUDE_DIRS "$ENV{NETCDF_ROOT}/include;$ENV{PNETCDF_ROOT}/include;${TPL_HDF5_INCLUDE_DIRS}"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(Netcdf_LIBRARY_DIRS "$ENV{NETCDF_ROOT}/lib;$ENV{PNETCDF_ROOT}/lib;${HDF5_LIBRARY_DIRS}"
  CACHE PATH "Set in ATTBDevEnv.cmake")
SET(Netcdf_LIBRARY_NAMES "netcdf;pnetcdf;${HDF5_LIBRARY_NAMES}"
  CACHE STRING "Set in ATTBDevEnv.cmake")

#
# Test disables
#

# See Trilinos #202
SET(STKUnit_tests_util_parallel_UnitTest_MPI_4_DISABLE ON
  CACHE BOOL  "Set in ATTBDevEnv.cmake")
# See Trilinos #211
SET(TeuchosNumerics_BLAS_ROTG_test_DISABLE ON
  CACHE BOOL  "Set in ATTBDevEnv.cmake")
SET(TeuchosNumerics_BLAS_ROTG_test_MPI_1_DISABLE ON
  CACHE BOOL  "Set in ATTBDevEnv.cmake")
