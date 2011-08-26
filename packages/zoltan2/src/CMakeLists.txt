
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE ZOLTAN2)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)
ASSERT_DEFINED(PROJECT_BINARY_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#
# Include packages/src directory
#
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Zoltan2_AlltoAll.hpp
  Zoltan2_Exceptions.hpp
  Zoltan2_IdentifierMap.hpp
  Zoltan2_Parameters.hpp
  Zoltan2_DebugManager.hpp
  Zoltan2_IdentifierMap_decl.hpp
  Zoltan2_IdentifierTraits.hpp
  Zoltan2_Partitioner.hpp
  Zoltan2_Environment.hpp
  Zoltan2_IdentifierMap_def.hpp
  Zoltan2_Memory.hpp
  Zoltan2_Util.hpp
  )

INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/tpetra/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/tpetra/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/teuchos/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/teuchos/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/kokkos/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/kokkos/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/kokkos/LinAlg)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/kokkos/NodeAPI)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/epetra/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/epetra/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/epetraext/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/epetraext/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/ThreadPool/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/ThreadPool/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/packages/triutils/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/packages/triutils/src)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR}/preCopyrightTrilinos/zoltan2/src)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/preCopyrightTrilinos/zoltan2/src)

# Set sources

APPEND_SET(SOURCES
  Zoltan2_Version.cpp
  Zoltan2_Util.cpp
  Zoltan2_Environment.cpp
  Zoltan2_DebugManager.cpp
  )

PACKAGE_ADD_LIBRARY(
  zoltan2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
