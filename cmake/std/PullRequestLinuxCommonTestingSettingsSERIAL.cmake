# This file contains the options needed to both run the pull request testing
# for Trilinos for the Linux GCC pull request testing builds, and to reproduce
# the errors reported by those builds. Prior to using this this file, the
# appropriate set of SEMS modules must be loaded and accessible through the
# SEMS NFS mount. (See the sems/PullRequestGCC*TestingEnv.sh files.)

# Usage: cmake -C PullRequestLinuxCommonTestingSettings.cmake

# Misc options typically added by CI testing mode in TriBITS

# Use the below option only when submitting to the dashboard
#set (CTEST_USE_LAUNCHERS ON CACHE BOOL "Set by default for PR testing")

set (Trilinos_TEST_CATEGORIES BASIC CACHE STRING "Set by default for PR testing")

set (Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES         ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ALLOW_NO_PACKAGES                    ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_DISABLE_ENABLED_FORWARD_DEP_PACKAGES ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_SECONDARY_TESTED_CODE         ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_IGNORE_MISSING_EXTRA_REPOSITORIES    ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_TESTS                         ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_CONFIGURE_TIMING              ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES      ON CACHE BOOL "Set by default for PR testing")
set (Trilinos_CTEST_USE_NEW_AAO_FEATURES           ON CACHE BOOL "Set by default for PR testing")

# Options from cmake/std/MpiReleaseDebugSharedPtSettings.cmake

set (CMAKE_BUILD_TYPE RELEASE CACHE STRING "Set by default for PR testing")

set (BUILD_SHARED_LIBS ON CACHE BOOL "Set by default for PR testing")

set (TPL_ENABLE_MPI                         OFF CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_DEBUG                  ON  CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_DEBUG_SYMBOLS          ON  CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_EXPLICIT_INSTANTIATION ON  CACHE BOOL "Set by default for PR testing")
set (Trilinos_ENABLE_SECONDARY_TESTED_CODE  OFF CACHE BOOL "Set by default for PR testing")
set (Teuchos_ENABLE_DEFAULT_STACKTRACE      OFF CACHE BOOL "Set by default for PR testing")

# Options from cmake/std/BasicCiTestingSettings.cmake

set (Trilinos_TPL_SYSTEM_INCLUDE_DIRS TRUE CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_GLM          OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Matio        OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_X11          OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Pthread      ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_BLAS         ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_LAPACK       ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Boost        ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_BoostLib     ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_METIS        OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_ParMETIS     OFF CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Zlib         ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_HDF5         ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_Netcdf       ON  CACHE BOOL "Set by default for PR testing")
set (TPL_ENABLE_SuperLU      ON  CACHE BOOL "Set by default for PR testing")
set (Trilinos_TRACE_ADD_TEST ON  CACHE BOOL "Set by default for PR testing")

set (TPL_ENABLE_Scotch OFF CACHE BOOL "Set by default for PR testing")

# Options from SEMSDevEnv.cmake
SET(CMAKE_C_COMPILER "$ENV{CC}" CACHE FILEPATH "Set by default for PR testing")

SET(CMAKE_CXX_COMPILER "$ENV{CXX}" CACHE FILEPATH "Set by default for PR testing")

SET(CMAKE_Fortran_COMPILER "$ENV{F90}" CACHE FILEPATH "Set by default for PR testing")

# SET(MPI_BASE_DIR "$ENV{SEMS_OPENMPI_ROOT}" CACHE PATH "Set by default for PR testing")

#SET(Trilinos_EXTRA_LINK_FLAGS "-lgomp -lgfortran -lldl -ldl" CACHE STRING "Set by default for PR testing")

SET(Trilinos_ENABLE_PyTrilinos OFF CACHE BOOL "Set by default for PR testing")

# Options (still from SEMSDevEnv.cmake) specific to TPLs
# TPL enables handled previously

SET(Boost_INCLUDE_DIRS "$ENV{SEMS_BOOST_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(Boost_LIBRARY_DIRS "$ENV{SEMS_BOOST_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(BoostLib_INCLUDE_DIRS "$ENV{SEMS_BOOST_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(BoostLib_LIBRARY_DIRS "$ENV{SEMS_BOOST_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(METIS_INCLUDE_DIRS "$ENV{SEMS_METIS_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(METIS_LIBRARY_DIRS "$ENV{SEMS_METIS_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

# SET(ParMETIS_INCLUDE_DIRS "$ENV{SEMS_PARMETIS_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
# SET(ParMETIS_LIBRARY_DIRS "$ENV{SEMS_PARMETIS_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(Zlib_INCLUDE_DIRS "$ENV{SEMS_ZLIB_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(Zlib_LIBRARY_DIRS "$ENV{SEMS_ZLIB_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(HDF5_INCLUDE_DIRS "$ENV{SEMS_HDF5_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(HDF5_LIBRARY_DIRS "$ENV{SEMS_HDF5_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(Netcdf_INCLUDE_DIRS "$ENV{SEMS_NETCDF_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(Netcdf_LIBRARY_DIRS "$ENV{SEMS_NETCDF_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

SET(SuperLU_INCLUDE_DIRS "$ENV{SEMS_SUPERLU_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
SET(SuperLU_LIBRARY_DIRS "$ENV{SEMS_SUPERLU_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")

# set (TPL_Scotch_INCLUDE_DIRS "$ENV{SEMS_SCOTCH_INCLUDE_PATH}" CACHE PATH "Set by default for PR testing")
# set (Scotch_LIBRARY_DIRS "$ENV{SEMS_SCOTCH_LIBRARY_PATH}" CACHE PATH "Set by default for PR testing")



