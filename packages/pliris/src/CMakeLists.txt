
INCLUDE(PackageLibraryMacros)
INCLUDE(CheckIncludeFileCXX)

#
# A) Package-specific configuration options
#

CHECK_INCLUDE_FILE_CXX(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE_CXX(time.h HAVE_TIME_H)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

APPEND_SET(HEADERS
  Pliris_ConfigDefs.h
  )

#
# Util files
#

APPEND_SET(HEADERS
  PlirisVersion.h
  )

APPEND_SET(SOURCES
  PlirisVersion.cpp
  )

IF (TPL_ENABLE_MPI)
  #
  # Extra Core files
  #

  APPEND_SET(HEADERS
    block.h
    distribute.h
    exchange_pivots.h
    factor.h
    forwardall.h
    forward.h
    pcomm.h
    perm1.h
    permute.h
    permute_mat.h
    permute_rhs.h
    solve.h
    vars.h
    x_factor.h
    x_solve.h
    xlu_solve.h
    )

  APPEND_SET(SOURCES
    distribute.c
    exchange_pivots.c
    factor_new.c
    forwardall.c
    forward.c
    pcomm.c
    permute.c
    permute_mat.c
    permute_rhs.c
    solve.c
    x_factor.c
    xlu_solve_new.c
    x_solve.c
    )

  #
  # Extra Util files
  #

  APPEND_SET(HEADERS
    BLAS_prototypes.h
    cblassp.h
    dblassp.h
    defines.h
    extern_vars.h
    macros.h
    mytime.h
    sblassp.h
    zblassp.h
    )

  APPEND_SET(SOURCES
    clean_code.c
    time.c
    )

# this is conditional on epetra in the makefile.am, however,
# epetra is a required dependency of pliris so it doesn't make
# sence to have this be conditional on something we must have
# to compile at all

  APPEND_SET(HEADERS
    Pliris.cpp
    )

  APPEND_SET(SOURCES
    Pliris.h
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
  ADD_DEFINITIONS(-DZCPLX)
  PACKAGE_ADD_LIBRARY(
    zpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()

IF(Pliris_ENABLE_SCPLX)
  ADD_DEFINITIONS(-DSCPLX)
  PACKAGE_ADD_LIBRARY(
    cpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()

IF(Pliris_ENABLE_DREAL)
  ADD_DEFINITIONS(-DDREAL)
  PACKAGE_ADD_LIBRARY(
    dpliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()

IF(Pliris_ENABLE_SREAL)
  ADD_DEFINITIONS(-SDREAL)
  PACKAGE_ADD_LIBRARY(
    spliris
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()
