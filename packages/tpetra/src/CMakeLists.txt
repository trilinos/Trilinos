
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Tpetra_CombineMode.hpp
  Tpetra_ConfigDefs.hpp
  Tpetra_CrsMatrix.hpp
  Tpetra_DefaultPlatform.hpp
  Tpetra_Directory.hpp
  Tpetra_DirectoryDecl.hpp
  Tpetra_DistObject.hpp
  Tpetra_Distributor.hpp
  Tpetra_Export.hpp
  Tpetra_Import.hpp
  Tpetra_ImportExportData.hpp
  Tpetra_Map.hpp
  Tpetra_MapData.hpp
  Tpetra_MapDataDecl.hpp
  Tpetra_MapDecl.hpp
  Tpetra_MultiVector.hpp
  Tpetra_MultiVectorData.hpp
  Tpetra_MultiVectorDataDecl.hpp
  Tpetra_MultiVectorDecl.hpp
  Tpetra_Operator.hpp 
  Tpetra_Platform.hpp
  Tpetra_SerialPlatform.hpp
  Tpetra_Util.hpp
  Tpetra_Vector.hpp
  Tpetra_VectorData.hpp
  Tpetra_VectorDataDecl.hpp
  Tpetra_VectorDecl.hpp
  Tpetra_Version.hpp
  )

SET(SOURCES ${SOURCES}
  Tpetra_Version.cpp
  )

IF (Trilinos_ENABLE_MPI)
  SET(HEADERS ${HEADERS}
    Tpetra_MpiPlatform.hpp
   )
ENDIF()

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  tpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
