
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

SET(HEADERS    "")
SET(SOURCES    "")
SET(SWIG_OPTS  "")
SET(INTERFACES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )


IF(PYTRILINOS_ENABLE_TEUCHOS)
ENDIF()

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS

  )

APPEND_SET(SOURCES

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

#PACKAGE_ADD_LIBRARY(
#  kokkos
#  HEADERS ${HEADERS}
#  SOURCES ${SOURCES}
#  )
