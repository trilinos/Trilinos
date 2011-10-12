
#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(CURRENT_PACKAGE ZOLTAN2)
CONFIGURE_FILE("${TRILINOS_HOME_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)

ASSERT_DEFINED(${PACKAGE_NAME}_HAVE_MALLINFO)
SET(HAVE_MALLINFO ${${PACKAGE_NAME}_HAVE_MALLINFO})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_METIS)
SET(HAVE_METIS ${${PACKAGE_NAME}_ENABLE_METIS})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ParMETIS)
SET(HAVE_PARMETIS ${${PACKAGE_NAME}_ENABLE_ParMETIS})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_PaToH)
SET(HAVE_PATOH ${${PACKAGE_NAME}_ENABLE_PaToH})

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Scotch)
SET(HAVE_SCOTCH ${${PACKAGE_NAME}_ENABLE_Scotch})


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
      Zoltan2_DebugManager.hpp
      Zoltan2_Environment.hpp
      Zoltan2_EpetraCrsGraphInput.hpp
      Zoltan2_EpetraCrsMatrixInput.hpp
      Zoltan2_Exceptions.hpp
      Zoltan2_GraphInput.hpp
      Zoltan2_GraphModel.hpp
      Zoltan2_IdentifierMap.hpp
      Zoltan2_IdentifierTraits.hpp
      Zoltan2_InputAdapter.hpp
      Zoltan2_MatrixInput.hpp
      Zoltan2_Memory.hpp
      Zoltan2_Model.hpp
      Zoltan2_Parameters.hpp
      Zoltan2_PartitioningProblem.hpp
      Zoltan2_PartitioningSolution.hpp
      Zoltan2_Problem.hpp
      Zoltan2_Solution.hpp
      Zoltan2_Standards.hpp
      Zoltan2_TemplateMacros.hpp
      Zoltan2_TpetraCrsGraphInput.hpp
      Zoltan2_TpetraCrsMatrixInput.hpp
      Zoltan2_Util.hpp
      Zoltan2_XpetraCrsGraphInput.hpp
      Zoltan2_XpetraCrsMatrixInput.hpp
  )

APPEND_SET(HEADERS
      algorithms/partition/Zoltan2_PartitioningAlgorithms.hpp
      algorithms/partition/Zoltan2_AlgScotch.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/partition)


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
  Zoltan2_Environment.cpp
  Zoltan2_Parameters.cpp
  Zoltan2_Version.cpp
  )


IF(${PACKAGE_NAME}_HAVE_MALLINFO)
APPEND_SET(SOURCES
  Zoltan2_Mallinfo.cpp
 )
ENDIF()

PACKAGE_ADD_LIBRARY(
  zoltan2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
