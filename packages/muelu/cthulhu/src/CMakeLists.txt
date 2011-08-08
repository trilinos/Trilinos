# CMakeLists inspired by Tpetra/src/CMakeLists.txt

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#
IF (${PACKAGE_NAME}_ENABLE_Tpetra)
  GLOBAL_SET(HAVE_CTHULHU_TPETRA ON)
  GLOBAL_SET(Cthulhu_ENABLE_Tpetra ON)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra)
  GLOBAL_SET(HAVE_CTHULHU_EPETRA ON)
  GLOBAL_SET(Cthulhu_ENABLE_Epetra ON)
ENDIF()

CONFIGURE_FILE(
  ../cmake/Cthulhu_config.hpp.in 
  ${CMAKE_CURRENT_BINARY_DIR}/Cthulhu_config.hpp
  )

# Note: cannot use the macro PACKAGE_CONFIGURE_FILE() because 
# Cthulhu_config.h.in is not in ${PACKAGE_SOURCE_DIR}/cmake
## PACKAGE_CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Cthulhu_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES 
  CrsMatrix/Cthulhu_Trans.cpp
  Comm/Cthulhu_Comm.cpp
  Map/Cthulhu_EpetraMap.cpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

# JG TODO: Is it the good macro to call ?
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Def)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Map)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Comm)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/VbrMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Vector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MultiVector)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/BlockMap)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Platform)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsGraph)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/DistObject)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Import)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/CrsMatrix)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Operator)

APPEND_GLOB(SOURCES */*.cpp)
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
  muelu-cthulhu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
