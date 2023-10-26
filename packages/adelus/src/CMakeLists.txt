
INCLUDE(CheckIncludeFileCXX)

#
# A) Package-specific configuration options
#

CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Core files
#

APPEND_SET(HEADERS
  Adelus_ConfigDefs.hpp
  )

#
# Util files
#

APPEND_SET(HEADERS
  AdelusVersion.hpp
  )

APPEND_SET(SOURCES
  AdelusVersion.cpp
  )

IF (TPL_ENABLE_MPI)
  #
  # Extra Core files
  #

  APPEND_SET(HEADERS
    Adelus_distribute.hpp
    Adelus_factor.hpp
    Adelus_forward.hpp
    Adelus_perm1.hpp
    Adelus_perm_mat.hpp
    Adelus_perm_rhs.hpp
    Adelus_solve.hpp
    Adelus_vars.hpp
    Adelus_xlu_solve.hpp
    Adelus_x_factor.hpp
    Adelus_x_solve.hpp
    Adelus_global_comm.hpp
    Adelus.hpp
    BlasWrapper_copy_spec.hpp
    BlasWrapper_copy.hpp
    )

  APPEND_SET(SOURCES
    Adelus_distribute.cpp
    )

  #
  # Extra Util files
  #

  APPEND_SET(HEADERS
    Adelus_defines.h
    Adelus_macros.h
    Adelus_mytime.hpp
    )

ENDIF()


#
# C) Define the targets for package's library(s)
#

# Using options to allow the compilation of all the possible libraries
# this is not what is done in the makefile.am, but what is done there
# is not really viable for cmake since all the libraries would be put on
# the link line for any package that depends on Adelus. That would cause
# pretty much every Adelus symbol to be multiply defined 4 times which
# isn't allowable. This seems like a reasonable comprimise

IF(Adelus_ENABLE_ZCPLX)
  MESSAGE("     Building the Z Complex library")
  IF(Adelus_ENABLE_TIMING)
    MESSAGE("     ... with timing")
  ENDIF()
  IF(Adelus_ENABLE_HOSTPINNEDMEM)
    MESSAGE("     ... use Cuda/HIP Host Pinned memory for MPI")
  ENDIF()
  IF(Adelus_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
  ENDIF()
  IF(Adelus_ENABLE_PRINTSTATUS)
    MESSAGE("     ... allow status print")
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    zadelus
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Adelus_ENABLE_SCPLX)
  MESSAGE("     Building the S Complex library")
  IF(Adelus_ENABLE_TIMING)
    MESSAGE("     ... with timing")
  ENDIF()
  IF(Adelus_ENABLE_HOSTPINNEDMEM)
    MESSAGE("     ... use Cuda/HIP Host Pinned memory for MPI")
  ENDIF()
  IF(Adelus_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
  ENDIF()
  IF(Adelus_ENABLE_PRINTSTATUS)
    MESSAGE("     ... allow status print")
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    cadelus
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Adelus_ENABLE_SREAL)
  MESSAGE("     Building the single precision library")
  IF(Adelus_ENABLE_TIMING)
    MESSAGE("     ... with timing")
  ENDIF()
  IF(Adelus_ENABLE_HOSTPINNEDMEM)
    MESSAGE("     ... use Cuda/HIP Host Pinned memory for MPI")
  ENDIF()
  IF(Adelus_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
  ENDIF()
  IF(Adelus_ENABLE_PRINTSTATUS)
    MESSAGE("     ... allow status print")
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    sadelus
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Adelus_ENABLE_DREAL)
  MESSAGE("     Building the double precision library")
  IF(Adelus_ENABLE_TIMING)
    MESSAGE("     ... with timing")
  ENDIF()
  IF(Adelus_ENABLE_HOSTPINNEDMEM)
    MESSAGE("     ... use Cuda/HIP Host Pinned memory for MPI")
  ENDIF()
  IF(Adelus_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
  ENDIF()
  IF(Adelus_ENABLE_PRINTSTATUS)
    MESSAGE("     ... allow status print")
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    dadelus
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ENDIF()
