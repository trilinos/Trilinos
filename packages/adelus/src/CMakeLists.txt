
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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Core files
#

#APPEND_SET(HEADERS
#  Pliris_ConfigDefs.h
#  )

#
# Util files
#

#APPEND_SET(HEADERS
#  PlirisVersion.h
#  )

#APPEND_SET(SOURCES
#  PlirisVersion.cpp
#  )

IF (TPL_ENABLE_MPI)
  #
  # Extra Core files
  #

  APPEND_SET(HEADERS
    block.h
    distribute.hpp
    factor.hpp
    pcomm.h
    perm1.hpp
    solve.hpp
    vars.hpp
    xlu_solve.hpp
    Pliris.hpp
    BlasWrapper_copy_spec.hpp
    BlasWrapper_copy.hpp
    BlasWrapper_iamax_spec.hpp
    BlasWrapper_iamax.hpp
    )

  APPEND_SET(SOURCES
    pcomm.c
    )

  #
  # Extra Util files
  #

  APPEND_SET(HEADERS
    defines.h
    macros.h
    mytime.hpp
    )

ENDIF()


#
# C) Define the targets for package's library(s)
#

# Using options to allow the compilation of all the possible libraries
# this is not what is done in the makefile.am, but what is done there
# is not really viable for cmake since all the libraries would be put on
# the link line for any package that depends on Pliris. That would cause
# pretty much every Pliris symbol to be multiply defined 4 times which
# isn't allowable. This seems like a reasonable comprimise

IF(Pliris_ENABLE_ZCPLX)
  MESSAGE("     Building the Z Complex library")
  ADD_DEFINITIONS(-DZCPLX)
  IF(Pliris_ENABLE_TIMING)
    MESSAGE("     ... with timing")
    ADD_DEFINITIONS(-DGET_TIMING)
  ENDIF()
  IF(Pliris_ENABLE_CUDAHOSTPINNED)
    MESSAGE("     ... use Cuda Host Pinned memory for MPI")
    ADD_DEFINITIONS(-DCUDA_HOST_PINNED_MPI)
  ENDIF()
  IF(Pliris_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
    ADD_DEFINITIONS(-DUSE_DEEPCOPY)
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    zpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Pliris_ENABLE_SCPLX)
  MESSAGE("     Building the S Complex library")
  ADD_DEFINITIONS(-DSCPLX)
  IF(Pliris_ENABLE_TIMING)
    MESSAGE("     ... with timing")
    ADD_DEFINITIONS(-DGET_TIMING)
  ENDIF()
  IF(Pliris_ENABLE_CUDAHOSTPINNED)
    MESSAGE("     ... use Cuda Host Pinned memory for MPI")
    ADD_DEFINITIONS(-DCUDA_HOST_PINNED_MPI)
  ENDIF()
  IF(Pliris_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
    ADD_DEFINITIONS(-DUSE_DEEPCOPY)
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    cpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Pliris_ENABLE_SREAL)
  MESSAGE("     Building the single precision library")
  ADD_DEFINITIONS(-DSREAL)
  IF(Pliris_ENABLE_TIMING)
    MESSAGE("     ... with timing")
    ADD_DEFINITIONS(-DGET_TIMING)
  ENDIF()
  IF(Pliris_ENABLE_CUDAHOSTPINNED)
    MESSAGE("     ... use Cuda Host Pinned memory for MPI")
    ADD_DEFINITIONS(-DCUDA_HOST_PINNED_MPI)
  ENDIF()
  IF(Pliris_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
    ADD_DEFINITIONS(-DUSE_DEEPCOPY)
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    spliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ELSEIF(Pliris_ENABLE_DREAL)
  MESSAGE("     Building the double precision(default) library")
  ADD_DEFINITIONS(-DDREAL)
  IF(Pliris_ENABLE_TIMING)
    MESSAGE("     ... with timing")
    ADD_DEFINITIONS(-DGET_TIMING)
  ENDIF()
  IF(Pliris_ENABLE_CUDAHOSTPINNED)
    MESSAGE("     ... use Cuda Host Pinned memory for MPI")
    ADD_DEFINITIONS(-DCUDA_HOST_PINNED_MPI)
  ENDIF()
  IF(Pliris_ENABLE_USEDEEPCOPY)
    MESSAGE("     ... use Kokkos::deep_copy for BLAS copy")
    ADD_DEFINITIONS(-DUSE_DEEPCOPY)
  ENDIF()
  TRIBITS_ADD_LIBRARY(
    dpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )

ENDIF()
