
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)
PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_InternalConfig.h)

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
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Stratimikos_DefaultLinearSolverBuilder.hpp
  Thyra_DefaultRealLinearSolverBuilder.hpp
  )

APPEND_SET(SOURCES
  Stratimikos_DefaultLinearSolverBuilder.cpp
  )

#
# C) Define the targets for package's library(s)
#

SET(DEPLIBS)

IF (${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(DEPLIBS stratimikosamesos)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)
  APPEND_SET(DEPLIBS stratimikosaztecoo)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(DEPLIBS stratimikosifpack)
ENDIF()

# ToDo: Add other solvers and preconditioners

PACKAGE_ADD_LIBRARY(
  stratimikos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )
