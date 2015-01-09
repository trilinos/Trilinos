# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(TrilinosCreateClientTemplateHeaders)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")

SET(SOURCES
  Utils/Xpetra_Utils.cpp
  Utils/Xpetra_EpetraUtils.cpp
  Platform/Xpetra_DefaultPlatform.cpp
  )
IF (${PACKAGE_NAME}_ENABLE_Epetra)
  APPEND_GLOB(SOURCES
    Map/Xpetra_EpetraMap.cpp
    Import/Xpetra_EpetraImport.cpp
    Export/Xpetra_EpetraExport.cpp
    MultiVector/Xpetra_EpetraMultiVector.cpp
    CrsGraph/Xpetra_EpetraCrsGraph.cpp
    CrsMatrix/Xpetra_EpetraCrsMatrix.cpp
    Vector/Xpetra_EpetraVector.cpp
    Vector/Xpetra_EpetraIntVector.cpp
    )
ENDIF()

#SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# mfh 18 Dec 2014: The "real" Kokkos(Classic) is now the TpetraClassic
# subpackage of Tpetra.  Thus, we have to test whether Tpetra is
# enabled, not whether Kokkos is enabled.  Furthermore, it wasn't
# previously correct to test if Kokkos was enabled, since Kokkos could
# have been enabled without the KokkosClassic subpackage having been
# enabled.
#
#IF (NOT ${PROJECT_NAME}_ENABLE_Kokkos)
IF (NOT ${PROJECT_NAME}_ENABLE_Tpetra)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/FakeKokkos)
  APPEND_GLOB(HEADERS FakeKokkos/*.hpp)
ENDIF()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/DistObject)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Export)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Import)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/RowMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils/ForwardDeclaration)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Vector)

#The HEADERS variable is the list of headers that will be copied to the install directory during "make install".
#Note: Do not use */*.hpp. Append FakeKokkos/*.hpp only if Kokkos is disabled.
APPEND_GLOB(HEADERS CrsGraph/*.hpp)
APPEND_GLOB(HEADERS CrsMatrix/*.hpp)
APPEND_GLOB(HEADERS DistObject/*.hpp)
APPEND_GLOB(HEADERS Export/*.hpp)
APPEND_GLOB(HEADERS Headers/*.hpp)
APPEND_GLOB(HEADERS Import/*.hpp)
APPEND_GLOB(HEADERS Map/*.hpp)
APPEND_GLOB(HEADERS MultiVector/*.hpp)
APPEND_GLOB(HEADERS Operator/*.hpp)
APPEND_GLOB(HEADERS Platform/*.hpp)
APPEND_GLOB(HEADERS RowGraph/*.hpp)
APPEND_GLOB(HEADERS RowMatrix/*.hpp)
APPEND_GLOB(HEADERS Utils/ForwardDeclaration/*.hpp)
APPEND_GLOB(HEADERS Utils/*.hpp)
APPEND_GLOB(HEADERS Vector/*.hpp)

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_SET(HEADERS ${DIR}/Xpetra_config.hpp)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
