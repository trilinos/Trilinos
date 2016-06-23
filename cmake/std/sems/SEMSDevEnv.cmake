#
# Base options for all SEMS Dev Env bulids for Trilinos
#

#
# A) Define the compilers
#
# NOTE: Set up different compilers depending on if MPI is enabled or not.
#

IF (TPL_ENABLE_MPI)
  # Set up MPI compiler wrappers
  ASSERT_DEFINED(ENV{MPICC})
  ASSERT_DEFINED(ENV{MPICXX})
  ASSERT_DEFINED(ENV{MPIF90})
  SET(CMAKE_C_COMPILER "$ENV{CC}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_Fortran_COMPILER "$ENV{FC}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
ELSE()
  MESSAGE(FATAL_ERROR  "ToDo: Set up support for serial compilers")
ENDIF()

# Add rpath for compiler libraries and gomp for parts built with OpenMP
IF (TPL_FIND_SHARED_LIBS)
  SET(LDL_LINK_ARG " -lldl")
ELSE()
  SET(LDL_LINK_ARG)
ENDIF()
SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS
  "-lgomp -lgfortran${LDL_LINK_ARG} -ldl"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# Point to the right MPI
ASSERT_DEFINED(ENV{SEMS_OPENMPI_ROOT})
SET(MPI_BASE_DIR "$ENV{SEMS_OPENMPI_ROOT}" CACHE PATH
 "Set in SEMSDevEnv.cmake")

# Turn on explicit template instantaition by default
SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION  ON  CACHE BOOL
  "Set in SEMSDevEnv.cmake")

# Turn off float and complex
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Sacado_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Thyra_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_DOUBLE OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_FLOAT OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Anasazi_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

# Enable configure timing
SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING ON CACHE BOOL
  "Set in SEMSDevEnv.cmake")

# Disable packages we can't enable on this system
SET(Trilinos_ENABLE_PyTrilinos OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

# Disable a bunch of TPLs that are not on this system
SET(TPL_ENABLE_GLM OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_Matio OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_SuperLU OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_X11 OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

# BLAS
SET(TPL_ENABLE_BLAS ON CACHE BOOL "Set in SEMSDevEnv.cmake")
# Assume BLAS is found in default path!

# LAPACK
SET(TPL_ENABLE_LAPACK ON CACHE BOOL "Set in SEMSDevEnv.cmake")
# Assume LAPACK is found in default path!

# Boost
SET(TPL_ENABLE_Boost ON CACHE BOOL "Set in SEMSDevEnv.cmake")
ASSERT_DEFINED(ENV{SEMS_BOOST_ROOT})
SET(Boost_INCLUDE_DIRS "$ENV{SEMS_BOOST_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Boost_LIBRARY_DIRS "$ENV{SEMS_BOOST_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")

# BoostLib
SET(TPL_ENABLE_BoostLib ON CACHE BOOL "Set in SEMSDevEnv.cmake")
ASSERT_DEFINED(ENV{SEMS_BOOST_ROOT})
SET(BoostLib_INCLUDE_DIRS "$ENV{SEMS_BOOST_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(BoostLib_LIBRARY_DIRS "$ENV{SEMS_BOOST_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")

# ZLIB
SET(TPL_ENABLE_ZLIB ON CACHE BOOL "Set in SEMSDevEnv.cmake")
ASSERT_DEFINED(ENV{SEMS_ZLIB_ROOT})
SET(TPL_ZLIB_INCLUDE_DIRS "$ENV{SEMS_ZLIB_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(ZLIB_LIBRARY_DIRS "$ENV{SEMS_ZLIB_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(ZLIB_LIBRARY_NAMES "z"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# HDF5
SET(TPL_ENABLE_HDF5 ON CACHE BOOL "Set in SEMSDevEnv.cmake")
ASSERT_DEFINED(ENV{SEMS_HDF5_ROOT})
SET(HDF5_INCLUDE_DIRS "$ENV{SEMS_HDF5_ROOT}/include;${TPL_ZLIB_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_DIRS "$ENV{SEMS_HDF5_ROOT}/lib;${ZLIB_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;${ZLIB_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# Netcdf
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "Set in SEMSDevEnv.cmake")
ASSERT_DEFINED(ENV{SEMS_NETCDF_ROOT})
SET(TPL_Netcdf_INCLUDE_DIRS "$ENV{SEMS_NETCDF_ROOT}/include;${TPL_HDF5_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Netcdf_LIBRARY_DIRS "$ENV{SEMS_NETCDF_ROOT}/lib;${HDF5_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Netcdf_LIBRARY_NAMES "netcdf;pnetcdf;${HDF5_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

#
# Test disables
#
