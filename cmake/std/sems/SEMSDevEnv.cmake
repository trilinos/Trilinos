#
# Base options for all SEMS Dev Env bulids for Trilinos
#


# Handle this being passed in with -C option instead of
# <Project>_CONFIGURE_OPTIONS_FILE.
IF ("${PROJECT_NAME}" STREQUAL "")
  SET(PROJECT_NAME Trilinos)
  FUNCTION(ASSERT_DEFINED VARS)
    FOREACH(VAR ${VARS})
      IF(NOT DEFINED ${VAR})
        MESSAGE(SEND_ERROR "Error, the variable ${VAR} is not defined!")
      ENDIF()
    ENDFOREACH()
  ENDFUNCTION()
ENDIF()

#
# A) Define the compilers and basic env
#

IF (TPL_ENABLE_MPI)
  # Set up MPI compiler wrappers
  ASSERT_DEFINED(ENV{MPICC})
  ASSERT_DEFINED(ENV{MPICXX})
  ASSERT_DEFINED(ENV{MPIF90})
  SET(CMAKE_C_COMPILER "$ENV{MPICC}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_CXX_COMPILER "$ENV{MPICXX}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_Fortran_COMPILER "$ENV{MPIF90}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  # Point to the right MPI
  ASSERT_DEFINED(ENV{SEMS_OPENMPI_ROOT})
  SET(MPI_BASE_DIR "$ENV{SEMS_OPENMPI_ROOT}" CACHE PATH
    "Set in SEMSDevEnv.cmake")
  # Make OpenMPI 1.6.5 mpi.h with GCC 4.8.3 warnings go away (see #1341)
  INCLUDE_DIRECTORIES(SYSTEM "$ENV{SEMS_OPENMPI_ROOT}/include")
  # NOTE: With TriBITS, all soruce files get built with the MPI compiler
  # wrappers so it should be safe to add this include directory with -isystem
  # to all builds.  The only worry is Fortran 90+ code that will cause
  # problems with older versions of gfortran.
ELSE()
  # Set up serial non-MPI compiler wrappers
  ASSERT_DEFINED(ENV{CC})
  ASSERT_DEFINED(ENV{CXX})
  ASSERT_DEFINED(ENV{F90})
  SET(CMAKE_C_COMPILER "$ENV{CC}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
  SET(CMAKE_Fortran_COMPILER "$ENV{F90}" CACHE FILEPATH
   "Set in SEMSDevEnv.cmake")
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

#
# B) Disable packages and TPLs by default not supported by SEMS Dev Env
#

# Don't have SWIG so can't enable PyTrilinos
SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF CACHE BOOL "Set in SEMSDevEnv.cmake")

#
# C) Set up the paths to the TPL includes and libs
#

SET(SEMS_MPI_VERSION $ENV{SEMS_MPI_VERSION})
#PRINT_VAR(SEMS_MPI_VERSION)

#
# Define helper function for finding the serial (non-MPI) version of a TPL for
# as serial build
#
# Called as:
#
#   SEMS_SELECT_TPL_ROOT_DIR( <semsTPLName> <tplRootDirOut>
#     [PARALLEL_EXT <parallelExt> SERIAL_EXT <serialExt>] )
#
# If arguments:
#
#   PARALLEL_EXT <parallelExt> SERIAL_EXT <serialExt>
#
# are not given, then it is assumed that <parallelExt>=parallel and
# <serialExt>=base.
#
FUNCTION(SEMS_SELECT_TPL_ROOT_DIR  SEMS_TPL_NAME  TPL_ROOT_DIR_OUT)

  CMAKE_PARSE_ARGUMENTS(
     #prefix
     PARSE
     # options
     ""
     #one_value_keywords
     ""
     #multi_value_keywords
     "PARALLEL_EXT;SERIAL_EXT"
     ${ARGN}
     )

  IF (PARSE_PARALLEL_EXT)
    SET(PARALLEL_EXT ${PARSE_PARALLEL_EXT})
    SET(SERIAL_EXT ${PARSE_SERIAL_EXT})
  ELSE()
    SET(PARALLEL_EXT "parallel")
    SET(SERIAL_EXT "base")
  ENDIF()

  SET(OPENMPI_VERSION_DIR "/openmpi/${SEMS_MPI_VERSION}/${PARALLEL_EXT}")
  #PRINT_VAR(OPENMPI_VERSION_DIR)
  
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
       #MESSAGE("-- " "This TPL is not pointing to a parallel version")
      SET(TPL_ROOT_DIR "${SEMS_TPL_ROOT_ENV_VAR}")
    ELSE()
      #MESSAGE("-- " "This TPL is pointing to a parallel version")
      STRING(REPLACE "${OPENMPI_VERSION_DIR}" "" TPL_ROOT_DIR_BASE
        "${SEMS_TPL_ROOT_ENV_VAR}")
      SET(TPL_ROOT_DIR "${TPL_ROOT_DIR_BASE}/${SERIAL_EXT}")
    ENDIF()
  ENDIF()

  #PRINT_VAR(TPL_ROOT_DIR)
  SET(${TPL_ROOT_DIR_OUT} "${TPL_ROOT_DIR}" PARENT_SCOPE)

ENDFUNCTION()

# Assume BLAS is found in default path!

# Assume LAPACK is found in default path!

# Boost
SEMS_SELECT_TPL_ROOT_DIR(BOOST Boost_ROOT)
#PRINT_VAR(Boost_ROOT)
SET(Boost_INCLUDE_DIRS "${Boost_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Boost_LIBRARY_DIRS "${Boost_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")

# BoostLib
SET(BoostLib_INCLUDE_DIRS "${Boost_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(BoostLib_LIBRARY_DIRS "${Boost_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")

# Scotch (SEMS only provides an MPI version)
IF (TPL_ENABLE_MPI)
  SEMS_SELECT_TPL_ROOT_DIR(SCOTCH Scotch_ROOT)
  SET(TPL_Scotch_INCLUDE_DIRS "${Scotch_ROOT}/include"
    CACHE PATH "Set in SEMSDevEnv.cmake")
  SET(Scotch_LIBRARY_DIRS "${Scotch_ROOT}/lib"
    CACHE PATH "Set in SEMSDevEnv.cmake")
ENDIF()

# ParMETIS (SEMS only provides an MPI version)
IF (TPL_ENABLE_MPI)
  SEMS_SELECT_TPL_ROOT_DIR(PARMETIS ParMETIS_ROOT)
  #PRINT_VAR(ParMETIS_ROOT)
  SET(TPL_ParMETIS_INCLUDE_DIRS "${ParMETIS_ROOT}/include"
    CACHE PATH "Set in SEMSDevEnv.cmake")
  SET(ParMETIS_LIBRARY_DIRS "${ParMETIS_ROOT}/lib"
    CACHE PATH "Set in SEMSDevEnv.cmake")
ENDIF()

# Zlib
SEMS_SELECT_TPL_ROOT_DIR(ZLIB Zlib_ROOT)
#PRINT_VAR(Zlib_ROOT)
SET(TPL_Zlib_INCLUDE_DIRS "${Zlib_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Zlib_LIBRARY_DIRS "${Zlib_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Zlib_LIBRARY_NAMES "z"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# HDF5
SEMS_SELECT_TPL_ROOT_DIR(HDF5 HDF5_ROOT)
#PRINT_VAR(HDF5_ROOT)
SET(HDF5_INCLUDE_DIRS "${HDF5_ROOT}/include;${TPL_Zlib_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_DIRS "${HDF5_ROOT}/lib;${Zlib_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(HDF5_LIBRARY_NAMES "hdf5_hl;hdf5;${Zlib_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# Netcdf
SEMS_SELECT_TPL_ROOT_DIR(NETCDF Netcdf_ROOT
  PARALLEL_EXT "exo_parallel" SERIAL_EXT "exo")
#PRINT_VAR(Netcdf_ROOT)
SET(TPL_Netcdf_INCLUDE_DIRS "${Netcdf_ROOT}/include;${TPL_HDF5_INCLUDE_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(Netcdf_LIBRARY_DIRS "${Netcdf_ROOT}/lib;${HDF5_LIBRARY_DIRS}"
  CACHE PATH "Set in SEMSDevEnv.cmake")
IF (TPL_ENABLE_MPI)
  SET(SEMS_PNETCDF_LIB_STR "pnetcdf;")
ELSE()
  SET(SEMS_PNETCDF_LIB_STR "")
ENDIF()
SET(Netcdf_LIBRARY_NAMES "netcdf;${SEMS_PNETCDF_LIB_STR};${HDF5_LIBRARY_NAMES}"
  CACHE STRING "Set in SEMSDevEnv.cmake")

# SuperLU
SEMS_SELECT_TPL_ROOT_DIR(SUPERLU SuperLU_ROOT)
#PRINT_VAR(SuperLU_ROOT)
SET(TPL_SuperLU_INCLUDE_DIRS "${SuperLU_ROOT}/include"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(SuperLU_LIBRARY_DIRS "${SuperLU_ROOT}/lib"
  CACHE PATH "Set in SEMSDevEnv.cmake")
SET(SuperLU_LIBRARY_NAMES "superlu;lapack;blas"
  CACHE STRING "Set in SEMSDevEnv.cmake")
