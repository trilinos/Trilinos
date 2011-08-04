
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

# We temporarily have a copy of MueLu::Cthulhu here.  Eventually it will be checked
# in and instead of having a copy we'll require it to build zoltan2.

SET(CTHULHU_DIR ${PROJECT_SOURCE_DIR}/preCopyrightTrilinos/zoltan2/cthulhu/src)

APPEND_SET(CTHULHU_SOURCES
   ${CTHULHU_DIR}/Comm/Cthulhu_Comm.cpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_Trans.cpp
   ${CTHULHU_DIR}/Def/Cthulhu_Debug.cpp
   ${CTHULHU_DIR}/Map/Cthulhu_EpetraMap.cpp
   ${CTHULHU_DIR}/Platform/Cthulhu_DefaultPlatform.cpp
   ${CTHULHU_DIR}/Vector/Cthulhu_CombineMode.cpp
  )

APPEND_SET(CTHULHU_HEADERS
   ${CTHULHU_DIR}/BlockMap/Cthulhu_BlockMap.hpp
   ${CTHULHU_DIR}/BlockMap/Cthulhu_EpetraBlockMap.hpp
   ${CTHULHU_DIR}/BlockMap/Cthulhu_TpetraBlockMap.hpp
   ${CTHULHU_DIR}/Comm/Cthulhu_Comm.hpp
   ${CTHULHU_DIR}/CrsGraph/Cthulhu_CrsGraph.hpp
   ${CTHULHU_DIR}/CrsGraph/Cthulhu_EpetraCrsGraph.hpp
   ${CTHULHU_DIR}/CrsGraph/Cthulhu_TpetraCrsGraph.hpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_CrsMatrix.hpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_CrsMatrixFactory.hpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_EpetraCrsMatrix.hpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_TpetraCrsMatrix.hpp
   ${CTHULHU_DIR}/CrsMatrix/Cthulhu_Trans.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_ConfigDefs.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_Debug.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_EpetraExceptions.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_Example.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_Exceptions.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_Parameters.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_UseDefaultTypes.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_UseShortNames.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_UseShortNamesOrdinal.hpp
   ${CTHULHU_DIR}/Def/Cthulhu_UseShortNamesScalar.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_EpetraExport.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_EpetraImport.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_Export.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_ExportFactory.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_Import.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_ImportFactory.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_TpetraExport.hpp
   ${CTHULHU_DIR}/Import/Cthulhu_TpetraImport.hpp
   ${CTHULHU_DIR}/Map/Cthulhu_EpetraMap.hpp
   ${CTHULHU_DIR}/Map/Cthulhu_Map.hpp
   ${CTHULHU_DIR}/Map/Cthulhu_MapFactory.hpp
   ${CTHULHU_DIR}/Map/Cthulhu_TpetraMap.hpp
   ${CTHULHU_DIR}/MultiVector/Cthulhu_EpetraMultiVector.hpp
   ${CTHULHU_DIR}/MultiVector/Cthulhu_MultiVector.hpp
   ${CTHULHU_DIR}/MultiVector/Cthulhu_MultiVectorFactory.hpp
   ${CTHULHU_DIR}/MultiVector/Cthulhu_TpetraMultiVector.hpp
   ${CTHULHU_DIR}/Operator/Cthulhu_CrsOperator.hpp
   ${CTHULHU_DIR}/Operator/Cthulhu_Operator.hpp
   ${CTHULHU_DIR}/Operator/Cthulhu_OperatorFactory.hpp
   ${CTHULHU_DIR}/Operator/Cthulhu_OperatorView.hpp
   ${CTHULHU_DIR}/Platform/Cthulhu_DefaultPlatform.hpp
   ${CTHULHU_DIR}/Platform/Cthulhu_MpiPlatform.hpp
   ${CTHULHU_DIR}/Platform/Cthulhu_SerialPlatform.hpp
   ${CTHULHU_DIR}/VbrMatrix/Cthulhu_TpetraVbrMatrix.hpp
   ${CTHULHU_DIR}/VbrMatrix/Cthulhu_VbrMatrix.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_CombineMode.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_EpetraIntVector.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_EpetraVector.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_TpetraVector.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_Vector.hpp
   ${CTHULHU_DIR}/Vector/Cthulhu_VectorFactory.hpp
  )


INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Comm)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/CrsMatrix)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Map)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Def)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/BlockMap)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/CrsGraph)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Import)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/MultiVector)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Operator)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Platform)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/VbrMatrix)
INCLUDE_DIRECTORIES(${CTHULHU_DIR}/Vector)



PACKAGE_ADD_LIBRARY(
  HEADERS ${HEADERS} ${CTHULHU_HEADERS}
  SOURCES ${SOURCES} ${CTHULHU_SOURCES}
  )
