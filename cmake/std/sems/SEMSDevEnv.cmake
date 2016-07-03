#
# Base options for all SEMS Dev Env bulids for Trilinos
#

#
# A) Set up some basic Trilinos options defaults
#

SET(BUILD_SHARED_LIBS ON CACHE BOOL
  "Set in SEMSDevEnv.cmake")

SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION ON CACHE BOOL
  "Set in SEMSDevEnv.cmake")

# Turn off float and complex by default
SET(Teuchos_ENABLE_FLOAT OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Teuchos_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Sacado_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Thyra_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_DOUBLE OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Tpetra_INST_COMPLEX_FLOAT OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(Anasazi_ENABLE_COMPLEX OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
# ToDo: Remove the above when Trlinos_ENABLE_FLOAT and Trilinos_ENABLE_COMPLEX
# are supported and are off by default (see Trilinos #362)

SET(${PROJECT_NAME}_ENABLE_CONFIGURE_TIMING ON CACHE BOOL
  "Set in SEMSDevEnv.cmake")

#
# B) Define the compilers and basic env
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
  ASSERT_DEFINED(ENV{SEMS_GCC_ROOT})
  SET(SEMS_GCC_ROOT "$ENV{SEMS_GCC_ROOT}")
  IF (SEMS_GCC_ROOT)
    SET(GCC_BIN "$ENV{SEMS_GCC_ROOT}/bin")
    #PRINT_VAR(GCC_BIN)
    IF ("${CMAKE_C_COMPILER}" STREQUAL "" OR "${CMAKE_C_COMPILER}" STREQUAL "gcc")
      SET(CMAKE_C_COMPILER "${GCC_BIN}/gcc" CACHE FILEPATH
       "Set in SEMSDevEnv.cmake" FORCE)
      # NOTE: For some reason, CMAKE_C_COMPILER can getset to "gcc" by default
      # which then finds the wrong gcc!
    ENDIF()
    SET(CMAKE_CXX_COMPILER "${GCC_BIN}/g++" CACHE FILEPATH
     "Set in SEMSDevEnv.cmake")
    SET(CMAKE_Fortran_COMPILER "${GCC_BIN}/gfortran" CACHE FILEPATH
     "Set in SEMSDevEnv.cmake")
  ELSE()
    MESSAGE(FATAL_ERROR
      "ERROR: Don't know how to set serial compilers for this env!")
  ENDIF()
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

#
# C) Disable packages and TPLs not supported by SEMS Dev Env
#

# Don't have SWIG so can't enable PyTrilinos
SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

# STK does not build in a serial build (see #466)
IF (NOT TPL_ENABLE_MPI)
  SET(${PROJECT_NAME}_ENABLE_STK OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
ENDIF()

# Disable a bunch of TPLs that are not on this system
SET(TPL_ENABLE_GLM OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_Matio OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_SuperLU OFF CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(TPL_ENABLE_X11 OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

# Disable Zoltan2 usage of 64-bit Scotch and ParMETIS becaues this causes
# several existing Zoltan2 tests to fail that pass otherwise (see Trilinos
# #476). Also, you have to not enable ParMETIS for ShyLU because it requires
# that ParMETIS is enabled in Zoltan2.
SET(Zoltan2_ENABLE_Scotch OFF CACHE BOOL "Disabled in SEAMSDevEnv.cmake")
SET(Zoltan2_ENABLE_ParMETIS OFF CACHE BOOL "Disabled in SEAMSDevEnv.cmake")
SET(ShyLUCore_ENABLE_ParMETIS OFF CACHE BOOL "Disabled in SEAMSDevEnv.cmake")

#
# D) Set up the paths to the TPL includes and libs
#

# Define helper function for finding the serial version of a TPL of this is a
# serial build.

SET(SEMS_MPI_VERSION $ENV{SEMS_MPI_VERSION})
#PRINT_VAR(SEMS_MPI_VERSION)
SET(OPENMPI_VERSION_DIR "/openmpi/${SEMS_MPI_VERSION}")
#PRINT_VAR(OPENMPI_VERSION_DIR)

FUNCTION(SEMS_SELECT_TPL_ROOT_DIR  SEMS_TPL_NAME  TPL_ROOT_DIR_OUT)
  #PRINT_VAR(SEMS_TPL_NAME)
  #PRINT_VAR(TPL_ROOT_DIR_OUT)
  SET(SEMS_TPL_ROOT_ENV_VAR_NAME SEMS_${SEMS_TPL_NAME}_ROOT)
  #PRINT_VAR(SEMS_TPL_ROOT_ENV_VAR_NAME)
  ASSERT_DEFINED(ENV{${SEMS_TPL_ROOT_ENV_VAR_NAME}})
  SET(SEMS_TPL_ROOT_ENV_VAR $ENV{${SEMS_TPL_ROOT_ENV_VAR_NAME}})
  #PRINT_VAR(SEMS_TPL_ROOT_ENV_VAR)
  IF (TPL_ENABLE_MPI)
    SET(TPL_ROOT_DIR "${SEMS_TPL_ROOT_ENV_VAR}")
  ELSE()
    # Serial build, so adjust the TPL root dir
    STRING(FIND "${SEMS_TPL_ROOT_ENV_VAR}" "${OPENMPI_VERSION_DIR}"
      OPENMPI_VERSION_DIR_IDX)
    #PRINT_VAR(OPENMPI_VERSION_DIR_IDX)
    IF ("${OPENMPI_VERSION_DIR_IDX}" STREQUAL "-1")
      # This TPL is not pointing to a parallel version
      SET(TPL_ROOT_DIR "${SEMS_TPL_ROOT_ENV_VAR}")
    ELSE()
      # This TPL is pointing to a parallel version
      STRING(REPLACE "${OPENMPI_VERSION_DIR}" "" TPL_ROOT_DIR_BASE
        "${SEMS_TPL_ROOT_ENV_VAR}")
      SET(TPL_ROOT_DIR "${TPL_ROOT_DIR_BASE}/base")
    ENDIF()
  ENDIF()
  SET(${TPL_ROOT_DIR_OUT} "${TPL_ROOT_DIR}" PARENT_SCOPE)
ENDFUNCTION()

# BLAS
SET(TPL_ENABLE_BLAS ON CACHE BOOL "Set in SEMSDevEnv.cmake")
# Above, assume BLAS is found in default path!

# LAPACK
SET(TPL_ENABLE_LAPACK ON CACHE BOOL "Set in SEMSDevEnv.cmake")
# Above, assume LAPACK is found in default path!

# Boost
SET(TPL_ENABLE_Boost ON CACHE BOOL "Set in SEMSDevEnv.cmake")
SEMS_SELECT_TPL_ROOT_DIR(BOOST Boost_ROOT)
#PRINT_VAR(Boost_ROOT)
SET(Boost_INCLUDE_DIRS "${Boost_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Boost_LIBRARY_DIRS "${Boost_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")

# BoostLib
SET(TPL_ENABLE_BoostLib ON CACHE BOOL "Set in SEMSDevEnv.cmake")
SET(BoostLib_INCLUDE_DIRS "${Boost_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(BoostLib_LIBRARY_DIRS "${Boost_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")


# Scotch (SEMS only provides an MPI version)
IF (TPL_ENABLE_MPI)
  # Disable 32-bit Scotch because it is not compatible with 64-bit ParMETIS
  # because it causes Zoltan Scotch tests to fail.
  #SET(TPL_ENABLE_Scotch ON CACHE BOOL "Set in SEMSDevEnv.cmake")
  #SEMS_SELECT_TPL_ROOT_DIR(SCOTCH Scotch_ROOT)
  #SET(TPL_Scotch_INCLUDE_DIRS "${Scotch_ROOT}/include"
  #  CACHE PATH "Set in SEMSDevEnv.cmake")
  #SET(Scotch_LIBRARY_DIRS "${Scotch_ROOT}/lib}"
  #  CACHE PATH "Set in SEMSDevEnv.cmake")
ENDIF()

# ParMETIS (SEMS only provides an MPI version)
IF (TPL_ENABLE_MPI)
  SET(TPL_ENABLE_ParMETIS ON CACHE BOOL "Set in SEMSDevEnv.cmake")
  SEMS_SELECT_TPL_ROOT_DIR(PARMETIS ParMETIS_ROOT)
  #PRINT_VAR(ParMETIS_ROOT)
  SET(TPL_ParMETIS_INCLUDE_DIRS "${ParMETIS_ROOT}/include"
    CACHE PATH "Set in SEMSDevEnv.cmake")
  SET(ParMETIS_LIBRARY_DIRS "${ParMETIS_ROOT}/lib}"
    CACHE PATH "Set in SEMSDevEnv.cmake")
ENDIF()

# Zlib
SET(TPL_ENABLE_Zlib ON CACHE BOOL "Set in SEMSDevEnv.cmake")
SEMS_SELECT_TPL_ROOT_DIR(ZLIB Zlib_ROOT)
#PRINT_VAR(Zlib_ROOT)
SET(TPL_Zlib_INCLUDE_DIRS "${Zlib_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Zlib_LIBRARY_DIRS "${Zlib_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Zlib_LIBRARY_NAMES "z"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# HDF5
SET(TPL_ENABLE_HDF5 ON CACHE BOOL "Set in SEMSDevEnv.cmake")
SEMS_SELECT_TPL_ROOT_DIR(HDF5 HDF5_ROOT)
#PRINT_VAR(HDF5_ROOT)
SET(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include;${TPL_Zlib_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_DIRS "${HDF5_ROOT}/lib;${Zlib_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;${Zlib_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# Netcdf
SET(TPL_ENABLE_Netcdf ON CACHE BOOL "Set in SEMSDevEnv.cmake")
SEMS_SELECT_TPL_ROOT_DIR(NETCDF Netcdf_ROOT)
#PRINT_VAR(Netcdf_ROOT)
SET(TPL_Netcdf_INCLUDE_DIRS "${Netcdf_ROOT}/include;${TPL_HDF5_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Netcdf_LIBRARY_DIRS "${Netcdf_ROOT}/lib;${HDF5_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Netcdf_LIBRARY_NAMES "netcdf;pnetcdf;${HDF5_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

#
# Test disables
#

# ToDo: Add test disables when needed!
