
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

APPEND_SET(HEADERS 
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Tpetra_ConfigDefs.hpp
  Tpetra_CrsGraph.hpp
  Tpetra_CrsMatrix.hpp
  Tpetra_CrsMatrixSolveOp.hpp
  Tpetra_DefaultPlatform.hpp
  Tpetra_Directory.hpp
  Tpetra_DirectoryDecl.hpp
  Tpetra_DistObject.hpp
  Tpetra_Distributor.hpp
  Tpetra_Export.hpp
  Tpetra_HybridPlatform.hpp
  Tpetra_Import.hpp
  Tpetra_ImportExportData.hpp
  Tpetra_Map.hpp
  Tpetra_MapDecl.hpp
  Tpetra_MultiVector.hpp
  Tpetra_MultiVectorDecl.hpp
  Tpetra_Operator.hpp 
  Tpetra_RowGraph.hpp
  Tpetra_RowMatrix.hpp
  Tpetra_SerialPlatform.hpp
  Tpetra_Util.hpp
  Tpetra_Vector.hpp
  Tpetra_VectorDecl.hpp
  Tpetra_Version.hpp
  )

# Set sources

APPEND_SET(SOURCES
  Tpetra_DefaultPlatform.cpp
  Tpetra_Distributor.cpp
  Tpetra_Version.cpp
  )

IF (TPL_ENABLE_MPI)
  SET(HEADERS ${HEADERS}
    Tpetra_MpiPlatform.hpp
   )
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
