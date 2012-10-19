
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Boost)
IF (${PACKAGE_NAME}_ENABLE_Boost)
  SET(HAVE_TEUCHOS_BOOST TRUE)
ENDIF()

ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}ParameterList)
IF (${PROJECT_NAME}_ENABLE_${PARENT_PACKAGE_NAME}ParameterList)
  SET(HAVE_TEUCHOS_PARAMETERLIST TRUE)
ENDIF()

TRIBITS_CONFIGURE_FILE(${PARENT_PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE TEUCHOS)
ASSERT_DEFINED(Trilinos_SOURCE_DIR)
CONFIGURE_FILE("${Trilinos_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h
ASSERT_DEFINED(Trilinos_BINARY_DIR)
INCLUDE_DIRECTORIES(${Trilinos_BINARY_DIR})
APPEND_SET(HEADERS
  ${Trilinos_BINARY_DIR}/Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_BINARY_DIR}/${PARENT_PACKAGE_NAME}_DLLExportMacro.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.h)
APPEND_GLOB(HEADERS ${DIR}/*.c)
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)

LIST(REMOVE_ITEM  HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Teuchos_SharedAccessPolicy.hpp
  )


#
# C) Define the targets for package's library(s)
#

IF (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  # Must suppress strict aliasing warnings for this file!
  SET_SOURCE_FILES_PROPERTIES(
    Teuchos_stacktrace.cpp
    PROPERTIES COMPILE_FLAGS -fno-strict-aliasing
    )
ENDIF()

TRIBITS_ADD_LIBRARY(
  teuchoscore
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEFINES -DTEUCHOS_LIB_EXPORTS_MODE
  )
