# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#
IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  GLOBAL_SET(HAVE_XPETRA_TPETRA ON)
  GLOBAL_SET(Xpetra_ENABLE_Tpetra ON)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  GLOBAL_SET(HAVE_XPETRA_EPETRA ON)
  GLOBAL_SET(Xpetra_ENABLE_Epetra ON)
ENDIF()

CONFIGURE_FILE(
  ../cmake/Xpetra_config.hpp.in 
  ${CMAKE_CURRENT_BINARY_DIR}/Xpetra_config.hpp
  )

# Note: cannot use the macro PACKAGE_CONFIGURE_FILE() because 
# Xpetra_config.h.in is not in ${PACKAGE_SOURCE_DIR}/cmake
## PACKAGE_CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Xpetra_config.h)

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

IF (NOT ${PACKAGE_NAME}_ENABLE_Kokkos)
  SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/FakeKokkos)
ENDIF()

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Def)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Utils)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Import)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Export)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/DistObject)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Vector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockedCrsOperator)

#APPEND_GLOB(SOURCES */*.cpp)
APPEND_GLOB(HEADERS */*.hpp)
#APPEND_GLOB(HEADERS ${DIR}/*/*.cpp)
#APPEND_GLOB(HEADERS ${DIR}/*/*.hpp)

PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
# APPEND_GLOB(HEADERS ${DIR}/*.hpp)
# APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu-xpetra
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
