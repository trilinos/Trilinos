

#
# A) Package-specific configuration options
#

IF (HAVE_STRATIMIKOS_AMESOS2)
  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  ASSERT_DEFINED(Amesos2_ENABLE_KLU2)
  IF (NOT (${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND Amesos2_ENABLE_KLU2))
    SET(HAVE_STRATIMIKOS_AMESOS2 OFF)
  ENDIF()
ENDIF()

# Set HAVE_BELOS_TPETRA if it's globally set by Belos package
# This is required for Stratimikos to register Tpetra-based Belos solvers
IF (DEFINED HAVE_BELOS_TPETRA)
  # Variable already set by Belos
ELSEIF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
  # Check if Belos has Tpetra support for backwards compatibility
  IF (TARGET Belos::belostpetra OR TARGET belostpetra)
    SET(HAVE_BELOS_TPETRA ON)
  ENDIF()
ENDIF()

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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Stratimikos_ConfigDefs.hpp
  Stratimikos_DefaultLinearSolverBuilder.hpp
  Stratimikos_LinearSolverBuilder_decl.hpp
  Stratimikos_LinearSolverBuilder_def.hpp
  Thyra_DefaultRealLinearSolverBuilder.hpp
  Stratimikos_InternalConfig.h
  )

APPEND_SET(SOURCES
  Stratimikos_LinearSolverBuilder.cpp
  )

INCLUDE(TrilinosCreateClientTemplateHeaders)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(HEADERS ${CMAKE_CURRENT_BINARY_DIR}/*.hpp)

#
# C) Define the targets for package's library(s)
#

SET(DEPLIBS)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos2)
IF (${PACKAGE_NAME}_ENABLE_ThyraTpetraAdapters AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  ASSERT_DEFINED(Amesos2_ENABLE_KLU2)
  IF (Amesos2_ENABLE_KLU2)
    APPEND_SET(DEPLIBS stratimikosamesos2)
  ENDIF()
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Belos)
IF (${PACKAGE_NAME}_ENABLE_Belos)
  APPEND_SET(DEPLIBS stratimikosbelos)
ENDIF()

# ToDo: Add other solvers and preconditioners

TRIBITS_ADD_LIBRARY(
  stratimikos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS ${DEPLIBS}
  )
