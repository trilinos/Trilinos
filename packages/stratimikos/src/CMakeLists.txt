
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_Config.h)
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_InternalConfig.h)

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
  Stratimikos_ConfigDefs.hpp
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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Amesos)
  APPEND_SET(DEPLIBS stratimikosamesos)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
IF (${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_AztecOO)
  APPEND_SET(DEPLIBS stratimikosaztecoo)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
IF (${PACKAGE_NAME}_ENABLE_Belos)
  APPEND_SET(DEPLIBS stratimikosbelos)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
IF (${PACKAGE_NAME}_ENABLE_Ifpack)
  APPEND_SET(DEPLIBS stratimikosifpack)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
IF (${PACKAGE_NAME}_ENABLE_ML)
  APPEND_SET(DEPLIBS stratimikosml)
ENDIF()

# ToDo: Add other solvers and preconditioners

TRIBITS_ADD_LIBRARY(
  stratimikos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )
