
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

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

#
# Core files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include)

APPEND_SET(HEADERS
  ../include/TPI.h
  ../include/TPI.hpp
  )

APPEND_SET(SOURCES
  TPI_pthreads.c
  TPI_Concurrency.c
  )

#
# Util files
#
APPEND_SET(SOURCES
  TPI_Walltime.c
  )

######################################

APPEND_SET(HEADERS
  )

APPEND_SET(SOURCES
  )

######################################
IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tpi
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
